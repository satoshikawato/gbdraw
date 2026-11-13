import { pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import { acceptNodeFeatureRecoveryRow, installSessionImportWorker } from './helpers/session-import-node.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { sessionLoadEvent, sessionFileOf, sessionLoadEventOfBytes, readSavedSession } from './helpers/session-documents.mjs';

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
installFakeSvgDom();
let restoredPrimaryReads = 0;
let restoredPrimaryArrayBufferReads = 0;
const restoredReadsByName = new Map();
const resourceReencodesByName = new Map();
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric({ name, resourceName }) {
    if (name === 'resourceReencodeCount') {
      resourceReencodesByName.set(
        resourceName,
        (resourceReencodesByName.get(resourceName) || 0) + 1
      );
    }
    // A restore reads a Session resource's bytes once (`readParsedSequenceRecords`); reading its text reads the
    // bytes first, so the byte read counts every content read.
    if (name !== 'resourceByteReadCount') return;
    restoredPrimaryReads += 1;
    restoredReadsByName.set(
      resourceName,
      (restoredReadsByName.get(resourceName) || 0) + 1
    );
  }
};
globalThis.File = class File extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = String(name || 'file');
    this.lastModified = options.lastModified ?? Date.now();
  }

  async text() {
    restoredPrimaryReads += 1;
    restoredReadsByName.set(
      this.name,
      (restoredReadsByName.get(this.name) || 0) + 1
    );
    return super.text();
  }

  async arrayBuffer() {
    restoredPrimaryArrayBufferReads += 1;
    return super.arrayBuffer();
  }
};
const alerts = [];
globalThis.alert = (message) => alerts.push(String(message));

installSessionImportWorker();
// A Session 27-44 moves to a Session 49 in Python, as the composition root's
// Worker port does it (DESIGN 4.3).
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const {
  adoptCanonicalRenderArtifacts,
  getCommittedCanonicalSession,
  serializeActiveRenderFiles,
  validateSessionLosatArtifacts,
  installShownDrawing
} = sessionConfig;
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { importSession, exportSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
const { promoteCanonicalRenderRequestToCurrent } = await import('../../gbdraw/web/js/services/session-request.js');
// The composition root's transform of an older Session's Results (R13 port).
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');

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

const circularSequenceCoverageSession = ({
  includeComparisonFasta = false,
  includeCatalogComparison = false
} = {}) => {
  const genbank = `LOCUS       REFERENCE                   8 bp    DNA     linear   UNK 01-JAN-1980
ACCESSION   REFERENCE
VERSION     REFERENCE.1
FEATURES             Location/Qualifiers
ORIGIN
        1 aaccggtt
//
`;
  const comparisonFasta = '>comparison\naaaacccc\n';
  /** @param {string} resourceId @param {string} name */
  const binding = (resourceId, name) => ({ resourceId, name, type: 'text/plain', lastModified: 0 });
  const renderRequest = {
    schema: 5,
    mode: 'circular',
    grouping: 'single',
    records: [{
      recordKey: 'reference-record',
      source: { kind: 'genbank', resourceId: 'reference-record' },
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
      pairwiseMatchStyle: 'ribbon',
      conservationBlastFiles: [{
        resourceId: 'comparison-blast',
        representation: 'file'
      }],
      ...(includeComparisonFasta
        ? {
            conservationFastaFiles: [{
              resourceId: 'comparison-fasta',
              representation: 'file'
            }]
          }
        : {}),
      conservationReference: 'subject',
      conservationLabels: ['Comparison'],
      conservationColors: ['#4e79a7'],
      conservationRingWidth: 14,
      conservationRingGap: 3
    },
    layout: {},
    comparisons: [],
    output: {
      prefix: 'circular-sequence-coverage',
      formats: ['interactive_svg'],
      overwrite: false,
      interactiveMetadataPolicy: 'auto'
    }
  };
  const referenceSource = {
    key: 'circular:record:0',
    recordId: 'REFERENCE.1',
    aliases: ['REFERENCE', 'REFERENCE.1'],
    sequence: 'AACCGGTT',
    origin: 'circular-reference',
    recordIndex: 0
  };
  const catalogComparisonSource = {
    key: 'homology:comparison:0:comparison',
    recordId: 'comparison',
    aliases: ['comparison'],
    sequence: 'AAAACCCC',
    origin: 'homology-comparison',
    sourceIndex: 0
  };
  return {
    format: 'gbdraw-session',
    version: 40,
    renderRequest,
    resources: {
      'reference-record': canonicalResource('genbank', 'reference.gbk', genbank),
      'comparison-blast': canonicalResource(
        'conservation-blast-file',
        'comparison.tsv',
        'comparison\tREFERENCE.1\t99\t4\t0\t0\t1\t4\t1\t4\t1e-20\t50\n'
      ),
      ...(includeComparisonFasta
        ? {
            'comparison-fasta': canonicalResource(
              'conservation-fasta-file',
              'comparison.fna',
              comparisonFasta
            )
          }
        : {})
    },
    // Every slot, as the Session 40 Web writer wrote its bindings
    // (`session-resources.js` at that version); a slot holds its file or null.
    webFiles: {
      bindings: {
        schema: 1,
        c_gb: binding('reference-record', 'reference.gbk'),
        c_gff: null,
        c_fasta: null,
        c_depth: null,
        c_conservation_blasts: [binding('comparison-blast', 'comparison.tsv')],
        c_conservation_blasts_source: null,
        c_conservation_fastas: [],
        c_conservation_sequence_sources: includeComparisonFasta ? [binding('comparison-fasta', 'comparison.fna')] : [],
        d_color: null,
        t_color: null,
        blacklist: null,
        whitelist: null,
        qualifier_priority: null,
        linearSeqs: [],
        linearComparisons: []
      }
    },
    config: {
      form: { legend: 'left' },
      adv: {
        circular_track_slots_enabled: false,
        circular_track_slots: [],
        circular_track_slots_axis_index: null,
        circular_track_slots_schema_version: 4,
        linear_track_slots_enabled: false,
        linear_track_slots: [],
        linear_track_slots_axis_index: null,
        linear_track_slots_schema_version: 2
      },
      circularConservation: {
        enabled: true,
        source: 'upload'
      }
    },
    ui: { mode: 'circular', selectedResultIndex: 0 },
    results: [{
      name: 'circular-sequence-coverage',
      content: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"><g id="diagram"></g></svg>'
    }],
    features: {},
    editorState: {
      featureCatalog: {
        schema: 3,
        items: [{
          resultIndex: 0,
          resultName: 'circular-sequence-coverage',
          recordKeys: ['reference-record'],
          features: [],
          biologicalFeatures: [],
          orthogroups: [],
          annotations: [],
          comparisonMatches: [{
            id: 'homology-comparison-1',
            match_kind: 'homology',
            query_record_id: 'comparison',
            subject_record_id: 'REFERENCE.1',
            source_index: 0,
            reference_side: 'subject',
            qstart: 1,
            qend: 4,
            sstart: 1,
            send: 4
          }],
          sequenceSources: [
            referenceSource,
            ...(includeCatalogComparison ? [catalogComparisonSource] : [])
          ]
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
};

const importSessionPayload = async (payload) => importSession({
  target: {
    files: [await sessionFileOf(payload)],
    value: 'selected'
  }
});

test('decision-required protein comparisons stay committed without projecting a partial draft', async () => {
  alerts.length = 0;
  // The resolved protein table's rows lie inside each record, as a writer
  // that resolved them against the records wrote them (PD-OI-073).
  const genbank = `LOCUS       TEST                       80 bp    DNA     linear   UNK 01-JAN-1980
FEATURES             Location/Qualifiers
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
       61 atgcatgcat gcatgcatgc
//
`;
  const proteinTable = `query\tsubject\tidentity\talignment_length\tmismatches\tgap_opens\tqstart\tqend\tsstart\tsend\tevalue\tbitscore
protein-a\tprotein-b\t95\t20\t1\t0\t10\t30\t50\t70\t1e-20\t120
`;
  // The typed resources as Python's encoder writes them (schema 1, the
  // schema of a Session 39 writer): `kind` names the value kind and the value
  // is a typed tree (`encode_canonical_typed_resource`).
  const emptyGroupFields = ['orthogroups', 'memberByProteinId', 'confidenceByOrthogroupId', 'descriptionsByOrthogroupId',
    'nameCandidatesByOrthogroupId', 'namesByOrthogroupId', 'orthologEdgesByOrthogroupId', 'orthologPathsByOrthogroupId',
    'rbhOrthogroups', 'relatedEdgesByOrthogroupId', 'scopeByOrthogroupId', 'sourceRecordIndexByOrthogroupId'];
  const orthogroupsJson = JSON.stringify({
    schema: 1,
    kind: 'orthogroupResult',
    value: { type: 'OrthogroupResult', fields: Object.fromEntries(emptyGroupFields.map((field) => [field, {}])) }
  });
  const collinearityJson = JSON.stringify({
    schema: 1,
    kind: 'blocks',
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
      ncbiBlastpBin: null,
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
  const file = await sessionFileOf(payload);
  const event = { target: { files: [file], value: 'selected' } };

  // Load promotes the held request with its context (C-OD14); the legacy
  // alignment target names a feature no catalog holds, so it is listed.
  const listed = [];
  const result = await importSession(event, {
    confirmLoadRepairs: async (dialog) => { listed.push(...dialog.items); return true; }
  });

  assert.equal(result.status, 'ok', JSON.stringify(result, (key, value) => (value instanceof Error ? value.message : value)));
  assert.deepEqual(result.repairs.map(({ code, path }) => ({ code, path })), [
    { code: 'request-part-dropped', path: 'renderRequest.layout.similarityAlignment' }
  ]);
  assert.equal(listed.length, 1);
  assert.equal(result.comparisonDisposition, 'DECISION_REQUIRED');
  assert.equal(state.activeDrawing().importedComparisonIntent.disposition, 'DECISION_REQUIRED');
  assert.equal(state.activeDrawing().importedComparisonIntent.action, null);
  assert.deepEqual(modeDrawing(state, 'linear').files.linearCanonicalComparisons, []);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges.length, 0);
  // The held comparisons are kept; the promotion moves only the alignment
  // target. Python's encoder writes the promoted request: it names the result
  // resources and writes the optional setting's default.
  const { alignOrthogroupFeature: _target, ...heldSettings } = generatedProteinComparison.settings;
  // Its decoder reads a released writer's "either" as "strand", the one it drew
  // as (OV-2301): the spelling of a retired value, not today's default.
  const promotedSettings = {
    ...heldSettings,
    collinearInferOrthogroups: true,
    collinearityParams: {
      ...heldSettings.collinearityParams,
      parameters: { ...heldSettings.collinearityParams.parameters, mergeOrientation: 'strand' }
    }
  };
  const committed = getCommittedCanonicalSession();
  assert.equal(committed.renderRequest.schema, CANONICAL_REQUEST_SCHEMA);
  const committedComparisons = committed.renderRequest.comparisons;
  assert.deepEqual(
    committedComparisons.map(({ resourceId: _id, ...comparison }) => comparison),
    [
      ...renderRequest.comparisons.slice(0, 3).map(({ resourceId: _id, ...comparison }) => comparison),
      { ...generatedProteinComparison, settings: promotedSettings }
    ]
  );
  // The table keeps its bytes; the typed results are in Python's current
  // schema, of the same kind.
  const committedData = committedComparisons.slice(0, 3).map(({ resourceId }) => committed.resources[resourceId].data);
  assert.equal(committedData[0], btoa(proteinTable));
  assert.deepEqual(
    committedData.slice(1).map((data) => JSON.parse(atob(data)).kind),
    ['orthogroupResult', 'blocks']
  );
  assert.equal(state.activeDrawing().losatProgram.value, 'blastn');
  assert.equal(state.losatExecution.executionMode, 'threaded');
  assert.equal(state.losatExecution.parallelWorkers, '3');
  assert.equal(state.losatExecution.totalThreadBudget, '12');
  assert.equal(state.activeDrawing().losat.blastn.task, 'dc-megablast');
  assert.equal(state.activeDrawing().losat.blastp.mode, 'collinear');
  assert.equal(state.activeDrawing().losat.blastp.maxHits, 1);
  assert.equal(state.selectedOrthogroupAlignmentFeature.value, '');

  modeDrawing(state, 'linear').files.linearCanonicalComparisons = [];
  // OV-3314: a candidate of the committed request (Rotate record to feature, Align, Align at feature) replaces only
  // the committed Session; the draft's comparisons stay as the draft has them, here the empty list a decision leaves.
  adoptCanonicalRenderArtifacts(getCommittedCanonicalSession(), {
    adoptOwnedRequest: true, drawing: modeDrawing(state, 'linear'), conservationPairs: 'committed'
  });
  assert.deepEqual(modeDrawing(state, 'linear').files.linearCanonicalComparisons, []);
  assert.equal(getCommittedCanonicalSession().renderRequest.comparisons.length, 4, 'the committed request keeps its comparisons');
  // Generate adopts what Load committed: the request of schema 9.
  adoptCanonicalRenderArtifacts({ renderRequest: committed.renderRequest, resources: committed.resources }, { drawing: modeDrawing(state, 'linear') });
  const collinearityName = committed.resources[committedComparisons[2].resourceId].name;
  const collinearityReencodesBefore = resourceReencodesByName.get(collinearityName) || 0;
  assert.deepEqual(
    modeDrawing(state, 'linear').files.linearCanonicalComparisons.map((comparison) => comparison.kind),
    [
      'precomputedProteinComparison',
      'orthogroupResult',
      'collinearityResult',
      'generatedProteinComparison'
    ]
  );

  const serialized = await serializeActiveRenderFiles(state, state.activeDrawing());
  assert.equal(serialized.linearCanonicalComparisons[0].file.data, btoa(proteinTable));
  assert.equal(serialized.linearCanonicalComparisons[1].file.data, committedData[1]);
  assert.equal(serialized.linearCanonicalComparisons[2].file.data, committedData[2]);
  assert.equal(
    resourceReencodesByName.get(collinearityName) || 0,
    collinearityReencodesBefore
  );
  assert.equal(serialized.linearCanonicalComparisons[2].valueKind, 'blocks');
  assert.deepEqual(
    serialized.linearCanonicalComparisons[3],
    { ...generatedProteinComparison, settings: promotedSettings }
  );
});

// D5: the page keeps no match sequences, so Load reads no input file for
// them, whether or not the catalog names a source for every match.
test('a current session loads without reading its input files for match sequences', async () => {
  alerts.length = 0;
  restoredPrimaryReads = 0;
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
      content: '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"><g id="diagram"></g></svg>'
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
  const file = await sessionFileOf(payload);
  const event = { target: { files: [file], value: 'selected' } };

  const result = await importSession(event);

  assert.equal(result.status, 'ok');
  assert.equal(state.drawingOf('linear').linearSeqs[0].gb.name, 'catalog-primary.gbk');
  assert.equal(restoredPrimaryReads, 0);
  assert.equal(restoredPrimaryArrayBufferReads, 0);

  const secondGenbank = genbank
    .replaceAll('CATALOG', 'SECOND')
    .replace('aaccggtt', 'ttggccaa');
  const sparsePayload = structuredClone(payload);
  sparsePayload.renderRequest.records.push({
    ...structuredClone(renderRequest.records[0]),
    recordKey: 'second-record',
    source: { kind: 'genbank', resourceId: 'second-record' }
  });
  sparsePayload.resources['second-record'] = canonicalResource(
    'genbank',
    'second-primary.gbk',
    secondGenbank
  );
  sparsePayload.webFiles.bindings.linearSeqs.push({
    ...structuredClone(sparsePayload.webFiles.bindings.linearSeqs[0]),
    uid: 'second-record',
    gb: {
      resourceId: 'second-record',
      name: 'second-primary.gbk',
      type: 'application/genbank',
      lastModified: 8
    }
  });
  sparsePayload.editorState.featureCatalog.items[0].recordKeys.push(
    'second-record'
  );
  restoredPrimaryReads = 0;
  restoredPrimaryArrayBufferReads = 0;
  const sparseResult = await importSession({
    target: {
      files: [await sessionFileOf(sparsePayload)],
      value: 'selected'
    }
  });

  assert.equal(sparseResult.status, 'ok');
  assert.equal(restoredPrimaryReads, 0);
  assert.equal(restoredPrimaryArrayBufferReads, 0);

  const missingPayload = structuredClone(sparsePayload);
  missingPayload.editorState.featureCatalog.items[0].comparisonMatches.push({
    id: 'missing-source-match',
    match_kind: 'pairwise',
    query_record_id: 'CATALOG.1',
    query_record_index: 0,
    subject_record_id: 'SECOND.1',
    subject_record_index: 1,
    qstart: 1,
    qend: 4,
    sstart: 1,
    send: 4
  });
  restoredPrimaryReads = 0;
  const missingResult = await importSession({
    target: {
      files: [await sessionFileOf(missingPayload)],
      value: 'selected'
    }
  });

  assert.equal(missingResult.status, 'ok');
  assert.equal(restoredPrimaryReads, 0);
});

test('a current Circular session loads without reading its reference or comparison files', async () => {
  restoredReadsByName.clear();
  const optionalResult = await importSessionPayload(
    circularSequenceCoverageSession()
  );

  assert.equal(optionalResult.status, 'ok');
  assert.equal(restoredReadsByName.get('reference.gbk') || 0, 0);
  assert.equal(restoredReadsByName.get('comparison.fna') || 0, 0);

  restoredReadsByName.clear();
  const recoverableResult = await importSessionPayload(
    circularSequenceCoverageSession({ includeComparisonFasta: true })
  );

  assert.equal(recoverableResult.status, 'ok');
  assert.equal(restoredReadsByName.get('reference.gbk') || 0, 0);
  assert.equal(restoredReadsByName.get('comparison.fna') || 0, 0);

  restoredReadsByName.clear();
  const catalogResult = await importSessionPayload(
    circularSequenceCoverageSession({ includeCatalogComparison: true })
  );

  assert.equal(catalogResult.status, 'ok');
  assert.equal(restoredReadsByName.get('reference.gbk') || 0, 0);
  assert.equal(restoredReadsByName.get('comparison.fna') || 0, 0);
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
  state.activeDrawing().circularConservation.series.splice(
    0,
    state.activeDrawing().circularConservation.series.length,
    {
      sourceKey: 'ring-0',
      sourceIndex: 4,
      label: 'Old label',
      color: '#000000',
      losat_gencode: 11
    },
    {
      sourceKey: 'ring-1',
      sourceIndex: 9,
      label: 'Second old label',
      color: '#111111',
      losat_gencode: 4
    }
  );

  // OV-2818: Show Comparison Rings turned off after the Generate stays off.
  state.activeDrawing().circularConservation.enabled = false;
  adoptCanonicalRenderArtifacts({
    renderRequest: promoteCanonicalRenderRequestToCurrent(renderRequest),
    resources,
    webFiles: {
      conservationBlastSource: 'losat-cache',
      conservationLosatFastaSources: [null, 'fasta']
    }
  }, { drawing: modeDrawing(state, 'circular'), conservationPairs: new Set(), installRings: true });

  assert.equal(state.activeDrawing().circularConservation.enabled, false);
  assert.equal(modeDrawing(state, 'circular').files.c_conservation_blasts_source, 'losat-cache');
  assert.equal(modeDrawing(state, 'circular').files.c_conservation_blasts.length, 2);
  assert.equal(modeDrawing(state, 'circular').files.c_conservation_blasts[0].size, 0);
  assert.equal(modeDrawing(state, 'circular').files.c_conservation_fastas.length, 2);
  // The flat form of a schema-5 request reads as one ring genome of one file, or of none.
  assert.deepEqual(modeDrawing(state, 'circular').files.c_conservation_fastas[0], []);
  assert.equal(modeDrawing(state, 'circular').files.c_conservation_fastas[1][0].name, 'comparison.fna');
  assert.equal(state.activeDrawing().circularConservation.source, 'losat');
  assert.deepEqual(
    {
      sourceIndex: state.activeDrawing().circularConservation.series[0].sourceIndex,
      label: state.activeDrawing().circularConservation.series[0].label,
      color: state.activeDrawing().circularConservation.series[0].color,
      losatGencode: state.activeDrawing().circularConservation.series[0].losat_gencode
    },
    {
      sourceIndex: 0,
      label: 'Missing FASTA',
      color: '#e15759',
      // The request carries no ring's genetic code, and the adoption takes nothing of the draft's ring at the same
      // index: the rings' adoption transition gives each ring the code of the draft's ring it is (OV-3178,
      // `adoptCircularRingResult`, pinned in circular-rings-setting.test.mjs).
      losatGencode: 1
    }
  );
  assert.equal(state.activeDrawing().circularConservation.series[1].sourceIndex, 1);

  installShownDrawing('circular');
  state.activeDrawing().circularConservation.enabled = true;
  const serialized = await serializeActiveRenderFiles(state, state.activeDrawing());
  assert.equal(serialized.c_conservation_blasts[0].data, '');
  assert.deepEqual(serialized.c_conservation_fastas[0], []);
  assert.equal(serialized.c_conservation_fastas[1][0].data, btoa(fasta));
});

// OV-3019: each ring Generate builds its request from the views the last one
// bound, adopts a new resource table and binds new views of its comparison
// FASTA. They stand for the bytes of the File the user chose, which rest on
// disk, so the History file estimate counts none of them (OV-2750), however
// many ring Generates the History keeps; that a payload whose bytes the page
// holds counts once is pinned in history-file-budget.test.mjs (OV-3019). A
// committed resource still names the bound view as its own (OV-1555).
test('ring Generates over one comparison FASTA count none of its bytes in the History file estimate (OV-3019, OV-2750)', async () => {
  const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
  const { resourceIdOwnedBy } = await import('../../gbdraw/web/js/services/resource-payload-owner.js');
  const { bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const genbank = 'LOCUS       REFERENCE                   8 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaccggtt\n//\n';
  const fasta = new File(['>comparison\naaaacccc\n'], 'comparison.fna', { type: 'text/plain' });
  installShownDrawing('circular');
  const drawing = state.activeDrawing();
  const inputs = drawing.slices.inputs;
  drawing.circularSources.splice(0, drawing.circularSources.length,
    { uid: 'record-1', gb: state.inputFiles.registerFile(new File([genbank], 'reference.gb')), gff: null, fasta: null });
  drawing.circularRecordList.value = [{ uid: 'record-1', k: 1, identity: 'record-1', selector: '#1', record_id: 'REFERENCE',
    record_length: 8, detectedTopology: 'linear' }];
  Object.assign(drawing.circularConservation, { enabled: true, source: 'losat', reference: 'subject', labels: '' });
  drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length,
    { sourceIndex: 0, fileName: fasta.name, label: 'comparison', color: '#4e79a7', losat_gencode: 11 });
  // One ring whose comparison genome is one file.
  bindDrawingInputFileTree(state.inputFiles, inputs, 'c_conservation_fastas', [[fasta]]);
  bindDrawingInputFileTree(state.inputFiles, inputs, 'c_conservation_sequence_sources', []);
  const fileIds = state.inputFiles.collectIds(drawing.files);
  // A ring Generate (run-analysis.js): its request names the comparison FASTA for the search and for the ring.
  const ringGenerate = async () => {
    const filesData = await serializeActiveRenderFiles(state, drawing);
    const canonical = buildCanonicalRenderRequest({
      state, drawing, filesData,
      resolvedCircularConservation: [{ sourceIndex: 0, name: 'comparison.tsv', text: '', label: 'comparison',
        color: '#4e79a7', fastas: filesData.c_conservation_fastas[0] }]
    });
    // A Generate that searched the rings installs them (`installRings`, which the rings' adoption transition decides).
    assert.equal(adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true, drawing, conservationPairs: new Set(), installRings: true }), true);
    state.inputFiles.collectIds(drawing.files, fileIds);
    return canonical;
  };
  await ringGenerate();
  await ringGenerate();
  const { resources, renderRequest } = await ringGenerate();
  const [[ringFasta]] = renderRequest.diagramOptions.conservationFastaFiles;
  const [[view]] = drawing.files.c_conservation_fastas;
  assert.ok(fileIds.size > 3, 'every ring Generate bound new views of the comparison FASTA');
  assert.equal(state.inputFiles.estimateBytes([...fileIds].filter((fileId) => (
    state.inputFiles.fileOf({ fileId })?.name === fasta.name)), sessionConfig.inputFilePayload), 0);
  assert.deepEqual(Object.keys(resources).filter((id) => resources[id].kind === 'conservation-fasta-file'), [ringFasta.resourceId],
    'the request holds the comparison FASTA once');
  assert.equal(resourceIdOwnedBy(resources, state.inputFiles.fileOf(view), 'conservation-fasta-file'), ringFasta.resourceId);
  assert.equal(resourceIdOwnedBy(resources, fasta, 'conservation-fasta-file'), ringFasta.resourceId);
});

// The Generate path of the app (`adoptCircularRingResult`, app-setup.js): only a Generate that searched rings installs
// them, then the rings reconcile with the rings it searched; a candidate of the committed request leaves them.
const { createCircularTrackSlotEditor } = await import('../../gbdraw/web/js/app/circular-track-slots.js');
const { createCircularConservationSources, createCircularConservationComparisonFiles } = await import('../../gbdraw/web/js/app/circular-conservation.js');
const ringKeyHelpers = await import('../../gbdraw/web/js/services/circular-conservation-keys.js');
const ringStackEditor = createCircularTrackSlotEditor({ state });
const ringSources = createCircularConservationSources({
  state,
  writeLosatCacheInfo: () => {},
  syncCircularConservationSlots: ringStackEditor.syncCircularConservationSlots,
  rekeyCircularConservationRows: ringStackEditor.rekeyCircularConservationRows,
  circularConservationRowsLeaving: ringStackEditor.circularConservationRowsLeaving,
  askingStackPorts: () => ({ askChoice: async () => null, runUndoable: (_label, fn) => fn(), busy: () => false, rows: () => [], mode: 'circular' })
});
// The shown Circular drawing's Remove of a ring (the comparison panel's ×) and Add Seq.
const ringFiles = createCircularConservationComparisonFiles({
  state,
  shownAction: (mode, action) => action(modeDrawing(state, mode)),
  writeLosatCacheInfo: () => {},
  runUndoable: async (_label, fn) => fn(),
  setRingGenomeLists: ringSources.setRingGenomeLists,
  setCircularConservationSourceFiles: ringSources.setCircularConservationSourceFiles,
  reconcileCircularConservation: ringSources.reconcileCircularConservation,
  refreshCircularConservationSlots: ringSources.refreshCircularConservationSlots,
  askRingSourceRemoval: ringSources.askRingSourceRemoval
});
/**
 * The app's adoption of a candidate (`adoptCanonicalRenderArtifacts` port) with Generate's own rings (`ringAdoption`,
 * null for a candidate of the committed request).
 */
const adoptRingResult = (canonical, options, ringAdoption) => ringSources.adoptCircularRingResult(options.drawing,
  (installRings) => adoptCanonicalRenderArtifacts(canonical, { ...options, installRings }), ringAdoption);
const ringFasta = (name, sequence = 'aaaacccc') => new File([`>${name}\n${sequence}\n`], name, { type: 'text/plain', lastModified: 7 });
const blastTableOf = (file) => `${file.name.replace(/\.[^.]+$/, '')}.circular_conservation.blastn.tsv`;
/** Each ring's fields, by label. */
const ringFields = (drawing) => Object.fromEntries(drawing.circularConservation.series.map(
  ({ label, color, schemeIndex, losat_gencode: gencode }) => [label, [color, schemeIndex, gencode]]));
/**
 * A Circular drawing whose LOSAT rings compare `files`, each ring's series entry with the fields `fields` sets (by
 * file place), the entries the reconcile gives them (real keys).
 * @param {File[]} files @param {Record<string, unknown>[]} fields @param {string} [colorScheme]
 */
const ringDrawing = async (files, fields, colorScheme = 'ten') => {
  const { isManagedConservationSlot, reconcileConservationSeries } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const { bindDrawingInputFileList, bindDrawingInputFileTree, setConservationBlastsSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const genbank = 'LOCUS       REFERENCE                   8 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaccggtt\n//\n';
  installShownDrawing('circular');
  const drawing = state.activeDrawing();
  drawing.circularSources.splice(0, drawing.circularSources.length,
    { uid: 'record-1', gb: state.inputFiles.registerFile(new File([genbank], 'reference.gb')), gff: null, fasta: null });
  drawing.circularRecordList.value = [{ uid: 'record-1', k: 1, identity: 'record-1', selector: '#1', record_id: 'REFERENCE',
    record_length: 8, detectedTopology: 'linear' }];
  Object.assign(drawing.circularConservation, { enabled: true, source: 'losat', reference: 'subject', labels: '', colorScheme });
  // No ring result bound: an earlier test's derived tables would make the rings a replay of them.
  bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', []);
  setConservationBlastsSource(drawing.slices.inputs, null);
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', files.map((file) => [file]));
  bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_sequence_sources', []);
  // The custom stack off and no ring row: an earlier test's ring rows would name rings this drawing does not have.
  drawing.adv.circular_track_slots_enabled = false;
  drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length,
    ...drawing.adv.circular_track_slots.filter((slot) => !isManagedConservationSlot(slot)));
  drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length,
    ...reconcileConservationSeries({ sourceFiles: files, previousSeries: [], legacyLabels: [] }).series
      .map((entry, index) => ({ ...entry, ...fields[index] })));
  return drawing;
};
/**
 * A ring Generate of the drawing (run-analysis.js): the rings it searches (`searchedRingSources`), its request
 * (`losatRequestRings`; each ring's table holds its genome file's name), and the app's adoption with those rings. `beforeAdoption` runs while it runs; a change of the
 * rings there makes the draft's rings no longer those it searched (`builtFromDraft`).
 * @returns {Promise<{ canonical: Record<string, any>, rings: Record<string, any>[] }>}
 */
const ringGenerate = async (drawing, tableName = blastTableOf, beforeAdoption = () => {}) => {
  const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
  const { ringGenomeFileLists } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const rings = ringKeyHelpers.searchedRingSources(drawing);
  const keysOf = () => JSON.stringify(ringKeyHelpers.searchedRingSources(drawing).map((ring) => ring.sourceKey));
  const before = keysOf();
  const filesData = await serializeActiveRenderFiles(state, drawing);
  const canonical = buildCanonicalRenderRequest({
    state, drawing, filesData,
    resolvedCircularConservation: ringKeyHelpers.losatRequestRings(rings, rings.map((ring) => ({ name: tableName(ring.file), text: `${ring.fileName}\n` })),
      ringGenomeFileLists(filesData.c_conservation_fastas))
  });
  beforeAdoption();
  adoptRingResult(canonical, { adoptOwnedRequest: true, drawing, conservationPairs: new Set() },
    { searched: rings, builtFromDraft: before === keysOf() });
  assert.equal(drawing.files.c_conservation_blasts_source, 'losat-cache', 'the Generate adopted a LOSAT-cache result');
  return { canonical, rings };
};
/** A candidate of the committed request that runs no ring search: a record rotation or an alignment. */
const adoptCommittedCandidate = (drawing) => adoptRingResult(getCommittedCanonicalSession(), {
  adoptOwnedRequest: true, drawing, conservationPairs: 'committed' }, null);
/** Save Session, then Load the file it wrote, with the app's ring reconcile after a Load; the loaded Circular drawing. */
const savedAndLoaded = async () => {
  Object.assign(globalThis.document, { body: { appendChild: () => {} },
    createElement: () => ({ click: () => {}, addEventListener: () => {}, parentNode: null }) });
  const saved = await exportSession('rings', { askChoice: async () => 'download', readKeptArtifacts: () => [] });
  assert.equal(saved?.status, 'saved', JSON.stringify(saved));
  const loaded = await importSession({
    target: { files: [new File([await saved.blob.arrayBuffer()], 'rings.gbdraw')], value: 'rings' }
  }, { confirmLoadRepairs: async () => true, installLoadedArtifactSlots: () => {} });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  ringSources.reconcileLoadedCircularConservation(state.activeDrawing(), true);
  return state.activeDrawing();
};

// PD-OI-110: a LOSAT ring without a colour of its own (Same color) sends '' as the CLI reads it, not an invented
// colour (OV-3235). The request holds only the colour each ring drew, so the adoption gives each ring the colour fields
// of the ring Generate searched (its own colour or none, its scheme index), as it gives its genetic code; a ring moved
// while the Generate ran keeps them with its genome.
test('a ring Generate sends no colour for a Same color ring and keeps each ring\'s fields with its genome', async () => {
  const { moveConservationSeriesEntry } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const drawing = await ringDrawing([ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')],
    [{ label: 'a', color: '#123456', schemeIndex: 1, losat_gencode: 1 }, { label: 'b', color: null, schemeIndex: 0, losat_gencode: 4 }],
    'same');
  const keyOf = Object.fromEntries(drawing.circularConservation.series.map(({ label, sourceKey }) => [label, sourceKey]));
  // The second ring moves outward while the Generate runs.
  const { canonical, rings } = await ringGenerate(drawing, blastTableOf,
    () => assert.ok(moveConservationSeriesEntry(drawing.circularConservation.series, 1, -1)));
  // The rings as Generate read them (`searchedRingSources`): the drawn colours and the series entries' colour fields.
  assert.deepEqual(rings.map(({ color, ownColor, schemeIndex, sourceKey }) => [color, ownColor, schemeIndex, sourceKey]),
    [['#123456', '#123456', 1, keyOf.a], ['', null, 0, keyOf.b]]);
  assert.deepEqual(canonical.renderRequest.diagramOptions.conservationColors, ['#123456', ''],
    'a Same color ring sends no colour (OV-3235)');
  assert.deepEqual(ringFields(drawing), { a: ['#123456', 1, 1], b: [null, 0, 4] },
    'each ring keeps its own entry\'s colour, scheme index and genetic code');

  // Later, with no Generate, the second ring moves outward again; then a candidate of the committed request (a record
  // rotation, an alignment) is adopted, which leaves the draft's rings as they are, so each genome keeps its colour and
  // scheme index.
  assert.ok(moveConservationSeriesEntry(drawing.circularConservation.series, 1, -1));
  assert.deepEqual(drawing.circularConservation.series.map(({ label }) => label), ['b', 'a']);
  adoptCommittedCandidate(drawing);
  assert.deepEqual(ringFields(drawing), { a: ['#123456', 1, 1], b: [null, 0, 4] }, 'a committed candidate keeps each ring\'s fields');
});

// PD-OI-110: each ring's colour fields (its own colour or none, its scheme index) and genetic code go through Save and
// Load with the ring, and a candidate after Load keeps them: also when two rings share a generated BLAST table name.

test('a loaded ring request keeps each ring\'s fields through a candidate with no edit', async () => {
  const drawing = await ringDrawing([ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')], [{ color: '#ff0000', losat_gencode: 4 }, {}]);
  await ringGenerate(drawing);
  const loaded = await savedAndLoaded();
  const before = ringFields(loaded);
  assert.deepEqual(before, { a: ['#ff0000', 0, 4], b: [null, 1, 1] });
  adoptCommittedCandidate(loaded);
  assert.deepEqual(ringFields(loaded), before);
});

test('a loaded ring request keeps each ring\'s fields after a ring file is added without a Generate', async () => {
  const drawing = await ringDrawing([ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')], [{ color: '#ff0000', losat_gencode: 4 }, {}]);
  await ringGenerate(drawing);
  const loaded = await savedAndLoaded();
  // Add Seq: the derived tables go, and each ring keeps its place on its FASTA files.
  ringFiles.addCircularConservationComparisonFile({ target: { files: [ringFasta('c.fna', 'aaaatttt')], value: 'c.fna' } });
  assert.deepEqual(loaded.circularConservation.series.map(({ fileName }) => fileName), ['a.fna', 'b.fna', 'c.fna']);
  adoptCommittedCandidate(loaded);
  assert.deepEqual(ringFields(loaded), { a: ['#ff0000', 0, 4], b: [null, 1, 1], c: [null, 2, 1] }, 'the rings keep their fields');
});

test('loaded rings that share a BLAST table name keep their own fields after they swap places', async () => {
  const { moveConservationSeriesEntry } = await import('../../gbdraw/web/js/services/conservation-series.js');
  // a.fna and a.gbk give one table name; so do two copies of one file.
  for (const [files, labels] of [
    [[ringFasta('a.fna'), ringFasta('a.gbk', 'aaaagggg')], ['a', 'a.gbk']],
    [[ringFasta('a.fna'), ringFasta('a.fna')], ['first', 'second']]
  ]) {
    const drawing = await ringDrawing(files, [{ label: labels[0], color: '#ff0000', losat_gencode: 4 }, { label: labels[1] }]);
    await ringGenerate(drawing, () => 'a.circular_conservation.blastn.tsv');
    const loaded = await savedAndLoaded();
    assert.deepEqual(ringFields(loaded), { [labels[0]]: ['#ff0000', 0, 4], [labels[1]]: [null, 1, 1] });
    assert.ok(moveConservationSeriesEntry(loaded.circularConservation.series, 1, -1));
    adoptCommittedCandidate(loaded);
    assert.deepEqual(ringFields(loaded), { [labels[0]]: ['#ff0000', 0, 4], [labels[1]]: [null, 1, 1] },
      `${files.map(({ name }) => name).join(' and ')}: the hand-set colour stays with its ring`);
  }
});

// Two copies of one file (same name, size and date), the second moved up before the Generate, or not: each copy keeps
// its own colour and genetic code through a candidate, Save, Load and a candidate after Load.
for (const [title, moveBeforeGenerate] of [['moved before the Generate', true], ['not moved', false]]) {
  test(`copies of one file keep their fields through Save, Load and a candidate: ${title}`, async () => {
    const { moveConservationSeriesEntry } = await import('../../gbdraw/web/js/services/conservation-series.js');
    const drawing = await ringDrawing([ringFasta('a.fna'), ringFasta('a.fna')],
      [{ label: 'first', color: '#ff0000', losat_gencode: 4 }, { label: 'second' }]);
    if (moveBeforeGenerate) assert.ok(moveConservationSeriesEntry(drawing.circularConservation.series, 1, -1));
    await ringGenerate(drawing, () => 'a.circular_conservation.blastn.tsv');
    const expected = moveBeforeGenerate
      ? { second: [null, 1, 1], first: ['#ff0000', 0, 4] }
      : { first: ['#ff0000', 0, 4], second: [null, 1, 1] };
    assert.deepEqual(ringFields(drawing), expected);
    adoptCommittedCandidate(drawing);
    assert.deepEqual(ringFields(drawing), expected, 'a candidate before Save');
    const loaded = await savedAndLoaded();
    assert.deepEqual(ringFields(loaded), expected, 'Load');
    adoptCommittedCandidate(loaded);
    assert.deepEqual(ringFields(loaded), expected, 'a candidate after Load');
  });
}

// A Session this writer saves after copies of one file moved since the Generate loads with the keys it saved, and each
// copy keeps its own colour and genetic code, also through a candidate after Load.
test('a Session whose copies moved after the Generate loads with its keys and each copy\'s fields', async () => {
  const { moveConservationSeriesEntry } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const expected = { second: [null, 1, 1], first: ['#ff0000', 0, 4] };
  const drawing = await ringDrawing([ringFasta('a.fna'), ringFasta('a.fna')],
    [{ label: 'first', color: '#ff0000', losat_gencode: 4 }, { label: 'second' }]);
  assert.ok(moveConservationSeriesEntry(drawing.circularConservation.series, 1, -1));
  await ringGenerate(drawing, () => 'a.circular_conservation.blastn.tsv');
  assert.ok(moveConservationSeriesEntry(drawing.circularConservation.series, 1, -1));
  const keysSaved = drawing.circularConservation.series.map(({ label, sourceKey }) => [label, sourceKey]);
  const loaded = await savedAndLoaded();
  assert.deepEqual(ringFields(loaded), expected, 'Load');
  assert.deepEqual(loaded.circularConservation.series.map(({ label, sourceKey }) => [label, sourceKey]), keysSaved,
    'a Session this change writes loads with the keys it saved');
  adoptCommittedCandidate(loaded);
  assert.deepEqual(ringFields(loaded), expected, 'a candidate after Load');
});

// The v0.14.5 Web app's own Save (tests/fixtures/sessions/README.md): TLOSATX rings a.fna, a.fna, b.fna, the second
// copy moved outward before the Generate. The tag kept the copy numbers of the file order, so the outer ring (the
// second copy, green, genetic code 11) is saved as copy 1 and the red one (genetic code 4) as copy 0 (OV-3238). Load
// matches copies of one file by the place each ring was stored at (`conservationSourcesOfSeries`), not by those
// numbers, so each ring keeps its colour and genetic code; so does a candidate of the committed request, and a Save and
// Load of this writer.
test('a v0.14.5 Session with two copies of a ring file moved before the Generate keeps each ring\'s fields', async () => {
  const fixture = readFileSync(new URL('../fixtures/sessions/web-losat-ring-copies.v50.gbdraw-session.json.gz', import.meta.url));
  const loaded = await importSession({
    target: { files: [new File([fixture], 'web-losat-ring-copies.v50.gbdraw-session.json.gz')], value: 'v0145' }
  }, { confirmLoadRepairs: async () => true, installLoadedArtifactSlots: () => {} });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  ringSources.reconcileLoadedCircularConservation(state.activeDrawing(), true);
  const placed = (drawing) => drawing.circularConservation.series.map(
    ({ sourceIndex, color, losat_gencode: gencode }) => [sourceIndex, color, gencode]);
  const expected = [[0, '#00ff00', 11], [1, '#ff0000', 4], [2, '#77b26f', 1]];
  assert.deepEqual(placed(state.activeDrawing()), expected, 'Load keeps each copy at its place with its fields');
  adoptCommittedCandidate(state.activeDrawing());
  assert.deepEqual(placed(state.activeDrawing()), expected, 'a candidate after Load');
  const resaved = await savedAndLoaded();
  assert.deepEqual(placed(resaved), expected, 'Save and Load of this writer');
  adoptCommittedCandidate(resaved);
  assert.deepEqual(placed(resaved), expected, 'a candidate after the second Load');
});

// A ring removed without a Generate stays removed through a candidate of the committed request (a rotation, an
// alignment), which leaves the draft's rings as they are (R1 (c)); the kept ring keeps its own colour, scheme index and
// genetic code. Two copies of one file, and two different files.
test('a ring kept by a Remove without a Generate keeps its fields through a candidate', async () => {
  for (const files of [[ringFasta('a.fna'), ringFasta('a.fna')], [ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')]]) {
    const drawing = await ringDrawing(files,
      [{ label: 'first', color: '#ff0000', losat_gencode: 4 }, { label: 'second', color: '#00ff00', losat_gencode: 11 }]);
    await ringGenerate(drawing, () => 'a.circular_conservation.blastn.tsv');
    ringFiles.removeCircularConservationSource(0);
    const what = files.map(({ name }) => name).join(' and ');
    assert.deepEqual(ringFields(drawing), { second: ['#00ff00', 1, 11] }, `${what}: Remove`);
    adoptCommittedCandidate(drawing);
    assert.deepEqual(ringFields(drawing), { second: ['#00ff00', 1, 11] }, `${what}: the candidate puts no ring back`);
  }
});

// OV-3239: a Remove without a Generate changes only the draft. The committed request keeps its own rings, each
// with its label, colour and table, through Save and Load; the next candidate of it (a record rotation, an alignment:
// `projectCommittedEditorIntent`) draws them so, and leaves the kept ring's fields in the draft.
test('a Remove without a Generate leaves the committed rings through Save, Load and a candidate (OV-3239)', async () => {
  const { projectCommittedEditorIntent, readCanonicalResourceText } = await import('../../gbdraw/web/js/services/session-request.js');
  // A loaded Session's tables live in its input copies, not in the descriptor: read each through the one reader.
  const drawnRings = (canonical) => {
    const options = canonical.renderRequest.diagramOptions;
    return Promise.all(options.conservationBlastFiles.map(async ({ resourceId }, index) => [options.conservationLabels[index],
      options.conservationColors[index], await readCanonicalResourceText(canonical.resources, resourceId)]));
  };
  for (const files of [[ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')], [ringFasta('a.fna'), ringFasta('a.fna')]]) {
    const drawing = await ringDrawing(files,
      [{ label: 'first', color: '#ff0000', losat_gencode: 4 }, { label: 'second', color: '#00ff00', losat_gencode: 11 }]);
    await ringGenerate(drawing, files[1].name === 'a.fna' ? () => 'a.circular_conservation.blastn.tsv' : blastTableOf);
    const committed = [['first', '#ff0000', 'a.fna\n'], ['second', '#00ff00', `${files[1].name}\n`]];
    assert.deepEqual(await drawnRings(getCommittedCanonicalSession()), committed);
    ringFiles.removeCircularConservationSource(0);
    const loaded = await savedAndLoaded();
    const what = files.map(({ name }) => name).join(' and ');
    assert.deepEqual(await drawnRings(getCommittedCanonicalSession()), committed, `${what}: the loaded committed request`);
    assert.deepEqual(await drawnRings(projectCommittedEditorIntent({ committed: getCommittedCanonicalSession(), state, drawing: loaded })),
      committed, `${what}: the candidate after Load`);
    adoptCommittedCandidate(loaded);
    assert.deepEqual(ringFields(loaded), { second: ['#00ff00', 1, 11] }, `${what}: the draft keeps the kept ring`);
  }
});

// OV-2442: a managed ring row of the custom stack names its ring by its source file's key, which a ring Generate
// changes (the derived BLAST table's). Each row moves with its ring and keeps its id and edits when a ring moved up
// before the Generate, also for copies of one file.
test('a ring Generate keeps each managed ring row with its ring', async () => {
  const { moveConservationSeriesEntry } = await import('../../gbdraw/web/js/services/conservation-series.js');
  for (const [files, tableName, whileRunning] of [
    [[ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')], blastTableOf, false],
    [[ringFasta('a.fna'), ringFasta('a.fna')], () => 'a.circular_conservation.blastn.tsv', false]
  ]) {
    const drawing = await ringDrawing(files, [{ label: 'first', color: '#ff0000' }, { label: 'second' }]);
    drawing.adv.circular_track_slots_enabled = true;
    ringStackEditor.syncCircularConservationSlots(drawing);
    const ringRows = () => drawing.adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation');
    ringRows().forEach((slot) => { slot.params = { ...slot.params, hand_edit: slot.id }; });
    const labelOfRow = () => Object.fromEntries(ringRows().map((slot) => [slot.id, [slot.params.label, slot.params.hand_edit]]));
    const before = labelOfRow();
    assert.deepEqual(Object.values(before).map(([label]) => label).sort(), ['first', 'second']);
    const moveUp = () => assert.ok(moveConservationSeriesEntry(drawing.circularConservation.series, 1, -1));
    if (!whileRunning) moveUp();
    const rowKeysBefore = ringRows().map((slot) => slot.params.series_key);
    await ringGenerate(drawing, tableName, whileRunning ? moveUp : () => {});
    assert.notDeepEqual(ringRows().map((slot) => slot.params.series_key), rowKeysBefore, 'the rows name the derived tables');
    assert.deepEqual(labelOfRow(), before,
      `${files.map(({ name }) => name).join(' and ')}, moved ${whileRunning ? 'while it ran' : 'before'}: each row stays with its ring, with its edits`);
  }
});

const ringRowFiles = (drawing) => drawing.adv.circular_track_slots
  .filter((slot) => slot.renderer === 'sequence_conservation').map((slot) => slot.params.fileName);
/**
 * A Circular drawing whose LOSAT rings compare `names` (no ring result bound), its custom stack in use from the
 * preset, then `arrange`d with the stack editor, under the ring colour scheme `colorScheme`.
 * @param {string[]} names @param {(drawing: DrawingState) => void} [arrange] @param {string} [colorScheme]
 */
const ringStackDrawing = async (names, arrange = () => {}, colorScheme = 'ten') => {
  const drawing = await ringDrawing(names.map((name, index) => ringFasta(`${name}.fna`, `aaaa${['cccc', 'gggg', 'tttt', 'cccg'][index]}`)),
    names.map((label) => ({ label })), colorScheme);
  ringSources.reconcileCircularConservation(drawing);
  ringStackEditor.setCircularTrackSlotsEnabled(true);
  ringStackEditor.resetCircularTrackSlotsToPreset(drawing.form.track_type);
  arrange(drawing);
  return drawing;
};

// A Rotate or an alignment after a Remove and an Add Seq without a Generate adopts the committed request, which leaves
// the draft's rings as they are (R1 (c)): the rows keep their ids, edits and places, and each ring its colour fields.
test('a committed candidate after a Remove and an Add Seq leaves the ring rows and colours', async () => {
  const drawing = await ringStackDrawing(['ringa', 'ringb', 'ringc']);
  await ringGenerate(drawing);
  const ringRows = () => drawing.adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation');
  ringFiles.removeCircularConservationSource(0);
  ringFiles.addCircularConservationComparisonFile({ target: { files: [ringFasta('ringz.fna', 'aaaacccg')], value: 'ringz.fna' } });
  ringRows().forEach((slot) => { slot.params = { ...slot.params, hand_edit: slot.id }; });
  const rows = () => ringRows().map((slot) => [slot.params.label, slot.id, slot.params.hand_edit]);
  const before = { rows: rows(), fields: ringFields(drawing) };
  assert.deepEqual(before.fields, { ringb: [null, 1, 1], ringc: [null, 2, 1], ringz: [null, 0, 1] }, 'ringz takes the freed index');
  adoptCommittedCandidate(drawing);
  assert.deepEqual({ rows: rows(), fields: ringFields(drawing) }, before);
  ringStackEditor.setCircularTrackSlotsEnabled(false);
});

// A reconcile after the Remove (an input-file Replace of the ring tables) keeps the kept ring's key, colour and genetic
// code, also through the next candidate.
test('a ring kept by a Remove keeps its fields through a Replace of the ring tables and a candidate', async () => {
  const { normalizeFileList } = await import('../../gbdraw/web/js/services/conservation-series.js');
  for (const files of [[ringFasta('a.fna'), ringFasta('a.fna')], [ringFasta('a.fna'), ringFasta('b.fna', 'aaaagggg')]]) {
    const drawing = await ringDrawing(files,
      [{ label: 'first', color: '#ff0000', losat_gencode: 4 }, { label: 'second', color: '#00ff00', losat_gencode: 11 }]);
    await ringGenerate(drawing, files[1].name === 'a.fna' ? () => 'a.circular_conservation.blastn.tsv' : blastTableOf);
    ringFiles.removeCircularConservationSource(0);
    const kept = drawing.circularConservation.series.map(({ label, sourceKey }) => [label, sourceKey]);
    ringSources.replaceCircularConservationBlasts(drawing, normalizeFileList(drawing.files.c_conservation_blasts));
    const what = files.map(({ name }) => name).join(' and ');
    assert.deepEqual(drawing.circularConservation.series.map(({ label, sourceKey }) => [label, sourceKey]), kept,
      `${what}: the Replace keeps the kept entry's key`);
    adoptCommittedCandidate(drawing);
    assert.deepEqual(ringFields(drawing).second, ['#00ff00', 1, 11], `${what}: the kept ring keeps its colour and genetic code`);
  }
});

// A ring row the stack lacks goes after the last ring row: an Add Seq, and a ring the series moved first while the
// stack was off when Use custom stack turns on. A Replace of the middle ring's table (uploaded, or the result of a ring
// search) keeps the ring at its place, and its row with it (`keepRingsAtPlaces`).
test('a ring row the stack lacks goes after the last ring row; a replaced table keeps its ring\'s row', async () => {
  const { moveConservationSeriesEntry } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const moveConservationSeries = (series) => [3, 2, 1].every((index) => moveConservationSeriesEntry(series, index, -1));
  const added = await ringStackDrawing(['ringa', 'ringb', 'ringc']);
  ringFiles.addCircularConservationComparisonFile({ target: { files: [ringFasta('ringz.fna', 'aaaacccg')], value: 'ringz.fna' } });
  assert.ok(moveConservationSeries(added.circularConservation.series));
  ringSources.refreshCircularConservationSlots(added);
  assert.deepEqual(ringRowFiles(added), ['ringa.fna', 'ringb.fna', 'ringc.fna', 'ringz.fna'], 'Add Seq, the ring moved first');
  ringStackEditor.setCircularTrackSlotsEnabled(false);

  const turnedOn = await ringStackDrawing(['ringa', 'ringb', 'ringc']);
  ringStackEditor.setCircularTrackSlotsEnabled(false);
  ringFiles.addCircularConservationComparisonFile({ target: { files: [ringFasta('ringz.fna', 'aaaacccg')], value: 'ringz.fna' } });
  assert.ok(moveConservationSeries(turnedOn.circularConservation.series));
  ringStackEditor.setCircularTrackSlotsEnabled(true);
  assert.deepEqual(ringRowFiles(turnedOn), ['ringa.fna', 'ringb.fna', 'ringc.fna', 'ringz.fna'], 'Use custom stack on after the move');
  ringStackEditor.setCircularTrackSlotsEnabled(false);

  const { setConservationBlastsSource, bindDrawingInputFileList } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const table = (name) => new File([`q\t${name}\t90\t8\t0\t0\t1\t8\t1\t8\t1e-5\t16\n`], name, { type: 'text/plain', lastModified: 7 });
  const uploaded = await ringStackDrawing([]);
  uploaded.circularConservation.source = 'upload';
  bindDrawingInputFileList(state.inputFiles, uploaded.slices.inputs, 'c_conservation_blasts', [table('a.tsv'), table('b.tsv'), table('c.tsv')]);
  setConservationBlastsSource(uploaded.slices.inputs, null);
  ringSources.reconcileCircularConservation(uploaded);
  assert.deepEqual(ringRowFiles(uploaded), ['a.tsv', 'b.tsv', 'c.tsv']);
  const [a, , c] = uploaded.files.c_conservation_blasts;
  ringSources.replaceCircularConservationBlasts(uploaded, [a, table('b2.tsv'), c]);
  assert.deepEqual(ringRowFiles(uploaded), ['a.tsv', 'b2.tsv', 'c.tsv'], 'a Replace of the middle ring table');
  ringStackEditor.setCircularTrackSlotsEnabled(false);
  uploaded.circularConservation.source = 'losat';

  const replayed = await ringStackDrawing(['ringa', 'ringb', 'ringc']);
  await ringGenerate(replayed);
  const [ringA, , ringC] = replayed.files.c_conservation_blasts;
  ringSources.replaceCircularConservationBlasts(replayed, [ringA, table('ringb2.circular_conservation.blastn.tsv'), ringC]);
  assert.deepEqual(ringRowFiles(replayed), ['ringa', 'ringb2', 'ringc'].map((name) => `${name}.circular_conservation.blastn.tsv`),
    'a Replace of the middle ring\'s result table');
  ringStackEditor.setCircularTrackSlotsEnabled(false);
});

test('legacy standalone config import migrates retired values without a writer envelope', async () => {
  alerts.length = 0;
  const file = await sessionFileOf({
    form: {
      prefix: 'legacy-standalone',
      linear_track_layout: 'tuckin'
    },
    adv: {
      label_placement: 'on_feature',
      multi_record_size_mode: 'sqrt'
    }
  });
  const event = { target: { files: [file], value: 'selected' } };

  const result = await importSession(event);

  assert.equal(result.status, 'legacy');
  assert.equal(state.activeDrawing().form.prefix, 'legacy-standalone');
  assert.equal(state.activeDrawing().form.linear_track_layout, 'below');
  assert.equal(state.activeDrawing().adv.label_placement, 'above_feature');
  assert.equal(state.activeDrawing().adv.multi_record_size_mode, 'auto');
  // The legacy notice belongs to the app; the service shows no native alert.
  assert.deepEqual(alerts, []);
});

const frozenV39Session = () => JSON.parse(gunzipSync(readFileSync(new URL(
  '../fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz',
  import.meta.url
))));

test('current schema-8 alignment retains exact typed anchors and explicit historical receipt absence on import', async () => {
  const source = await readSavedSession(new Blob([readFileSync(new URL(
    '../test_inputs/BGC0000708-BGC0000713.gbdraw',
    import.meta.url
  ))]));
  const file = await sessionFileOf(source);
  const result = await importSession({ target: { files: [file], value: 'selected' } });

  assert.equal(result.status, 'ok');
  assert.equal(state.similarityAlignmentPlan.value.groupId, 'og_1');
  assert.deepEqual(state.similarityAlignmentPlan.value.reference, {
    recordKey: source.drawings.find((entry) => entry.id === source.ui.activeDrawingId).renderRequest.records[0].recordKey,
    biologicalFeatureId: 'fed46a3a6',
    sourceFeatureIndex: 14,
    stableFeatureSvgId: 'fed46a3a6'
  });
  assert.equal(state.linearRecordTranslations.value.length, 5);
});

test('visible protein cache entries from the frozen v39 session recover stable edge identities', async () => {
  alerts.length = 0;
  const sessionData = frozenV39Session();
  const sourcePairs = sessionData.losatCache.entries
    .filter((entry) => entry.display !== false)
    .map((entry) => [entry.queryRecordInstanceKey, entry.subjectRecordInstanceKey]);
  const file = await sessionFileOf(sessionData);
  const result = await importSession(
    { target: { files: [file], value: 'selected' } },
    { transformLegacyResultSvg, confirmLoadRepairs: acceptNodeFeatureRecoveryRow }
  );

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

// OV-987: a Collinear display search can run in the other orientation than the
// edge whose card shows it (a CLI grid that places a record above one listed
// before it). Save writes the search's record instance keys and the edge's
// fields; Load keeps the edge, as Generate does.
test('a visible cache entry stored against its edge orientation keeps the edge at Load', async () => {
  alerts.length = 0;
  const sessionData = frozenV39Session();
  const entries = sessionData.losatCache.entries;
  const [edgeEntry] = entries.filter((entry) => entry.display !== false);
  const reversed = entries.find((entry) => (
    entry.queryRecordInstanceKey === edgeEntry.subjectRecordInstanceKey
    && entry.subjectRecordInstanceKey === edgeEntry.queryRecordInstanceKey
  ));
  const edgeQueryUid = edgeEntry.queryRecordInstanceKey;
  const edgeSubjectUid = edgeEntry.subjectRecordInstanceKey;
  edgeEntry.display = false;
  Object.assign(reversed, {
    display: true,
    filename: 'reversed.tsv',
    edgeKey: `${edgeQueryUid}->${edgeSubjectUid}`,
    queryUid: edgeQueryUid,
    subjectUid: edgeSubjectUid,
    queryIndex: 0,
    subjectIndex: 1,
    ordinal: 0
  });
  const file = await sessionFileOf(sessionData);
  const result = await importSession(
    { target: { files: [file], value: 'selected' } },
    { transformLegacyResultSvg, confirmLoadRepairs: acceptNodeFeatureRecoveryRow }
  );

  assert.equal(result.status, 'ok');
  const row = state.losatCacheInfo.value.find((entry) => entry.key === reversed.key);
  assert.deepEqual(
    [row?.edgeKey, row?.queryUid, row?.subjectUid, row?.queryIndex, row?.subjectIndex, row?.ordinal],
    [`${edgeQueryUid}->${edgeSubjectUid}`, edgeQueryUid, edgeSubjectUid, 0, 1, 0]
  );
});

// OV-987, CLI-written entry: the CLI stores only the search's record instance
// keys; a display search in the other orientation than its edge takes the edge.
test('a CLI-written visible cache entry stored against its edge orientation takes the edge at Load', async () => {
  alerts.length = 0;
  const sessionData = frozenV39Session();
  const entries = sessionData.losatCache.entries;
  const [edgeEntry] = entries.filter((entry) => entry.display !== false);
  const reversed = entries.find((entry) => (
    entry.queryRecordInstanceKey === edgeEntry.subjectRecordInstanceKey
    && entry.subjectRecordInstanceKey === edgeEntry.queryRecordInstanceKey
  ));
  const edgeQueryUid = edgeEntry.queryRecordInstanceKey;
  const edgeSubjectUid = edgeEntry.subjectRecordInstanceKey;
  edgeEntry.display = false;
  reversed.display = true;
  for (const field of ['edgeKey', 'queryUid', 'subjectUid', 'queryIndex', 'subjectIndex', 'ordinal']) {
    assert.equal(reversed[field], undefined, field);
  }
  const file = await sessionFileOf(sessionData);
  const result = await importSession(
    { target: { files: [file], value: 'selected' } },
    { transformLegacyResultSvg, confirmLoadRepairs: acceptNodeFeatureRecoveryRow }
  );

  assert.equal(result.status, 'ok');
  const row = state.losatCacheInfo.value.find((entry) => entry.key === reversed.key);
  assert.deepEqual(
    [row?.edgeKey, row?.queryUid, row?.subjectUid, row?.queryIndex, row?.subjectIndex, row?.ordinal],
    [`${edgeQueryUid}->${edgeSubjectUid}`, edgeQueryUid, edgeSubjectUid, 0, 1, 0]
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
  const file = await sessionFileOf(sessionData);
  const result = await importSession(
    { target: { files: [file], value: 'selected' } },
    { transformLegacyResultSvg, confirmLoadRepairs: acceptNodeFeatureRecoveryRow }
  );

  assert.equal(result.status, 'ok');
  assert.equal(state.losatCacheInfo.value[0].edgeKey, undefined);
  assert.equal(
    state.losatCacheInfo.value.filter((entry) => entry.edgeKey === conflictingEdgeKey).length,
    1
  );
});

// A rejected Session artifact is a classified Session-schema diagnostic (R6).
const sessionSchemaDiagnostic = (error) => {
  assert.equal(error.code, 'INPUT_INVALID');
  assert.deepEqual(error.context, { field: 'schema', reason: 'FIELDS' });
  return true;
};

const syntheticConservationSession = () => JSON.parse(gunzipSync(readFileSync(new URL(
  '../fixtures/sessions/synthetic_conservation.gbdraw-session.json.gz',
  import.meta.url
))).toString('utf8'));

// PD-OI-110 Must preserve: a Session up to 50 stored a colour for every ring (its Web writer filled them), so each
// is the ring's own. The rings open, show in the panel, and Generate draws them in the colours they were saved with.
test('a Session 40 with ring colours opens, shows and generates them unchanged', async () => {
  const { effectiveRingColor, orderedConservationSources } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
  const stored = ['#4e79a7', '#e15759', '#59a14f'];
  const loaded = await importSession({
    target: { files: [new Blob([JSON.stringify(syntheticConservationSession())], { type: 'application/json' })], value: 'selected' }
  });
  assert.equal(loaded.status, 'ok');
  const drawing = modeDrawing(state, 'circular');
  const { circularConservation } = drawing;
  assert.equal(circularConservation.colorScheme, 'same');
  assert.deepEqual(circularConservation.series.map((entry) => entry.color), stored, 'each colour is the ring\'s own');
  assert.deepEqual(circularConservation.series.map((entry) => effectiveRingColor(entry, circularConservation.colorScheme)), stored,
    'the panel swatches');
  assert.deepEqual(orderedConservationSources(drawing.files.c_conservation_blasts, circularConservation).map(({ color }) => color),
    stored, 'the rings Generate reads');
  installShownDrawing('circular');
  const filesData = await serializeActiveRenderFiles(state, drawing);
  assert.deepEqual(buildCanonicalRenderRequest({ state, drawing, filesData }).renderRequest.diagramOptions.conservationColors, stored,
    'the request Generate sends');
});

// E-OD02: Load drops the repeat of a cache key with an unlisted row
// (E-OD05) and loads the rest; Save refuses the repeat
// (`validateSessionLosatArtifacts`, the tests below).
test('a Session that repeats a raw LOSAT cache key loads without the repeat', async () => {
  const valid = syntheticConservationSession();
  const accepted = await importSession(await sessionLoadEvent(valid));
  assert.equal(accepted.status, 'ok');
  const validKeys = state.losatCacheInfo.value.map((entry) => entry.key);

  alerts.length = 0;
  const repeated = syntheticConservationSession();
  const [, ringC] = repeated.losatCache.entries;
  repeated.losatCache.entries.push({
    ...ringC,
    filename: 'comparison-c-copy.circular_conservation.losatn.tsv'
  });
  const result = await importSession(await sessionLoadEvent(repeated));
  assert.equal(result.status, 'ok');
  assert.deepEqual(result.repairs.map((repair) => repair.code), ['cache-entry-dropped']);
  // The first entry under the key stays; the repeat is dropped.
  assert.deepEqual(state.losatCacheInfo.value.map((entry) => entry.key), validKeys);
  assert.deepEqual(alerts, []);
});

// OV-981: a saved raw LOSAT file name that is no text is read as a missing
// one (the entry's default name), so a download is never named
// "[object Object]". The entry and its rows load; no row is listed.
test('a Session raw LOSAT cache file name that is no text loads as the default name', async () => {
  const session = syntheticConservationSession();
  const [first, second] = session.losatCache.entries;
  first.filename = {};
  second.filename = 7;
  alerts.length = 0;
  const result = await importSession(await sessionLoadEvent(session));
  assert.equal(result.status, 'ok');
  assert.deepEqual(result.repairs, []);
  assert.deepEqual(state.losatCacheInfo.value.map((entry) => entry.filename), [
    'losat_pair_1.tsv', 'losat_pair_2.tsv', 'comparison-d.circular_conservation.losatn.tsv'
  ]);
  assert.deepEqual(alerts, []);
});

// OV-983: a search context that is no text (a list of the 64 hex characters)
// is no current raw entry: Load drops it with an unlisted row, as any entry
// it cannot read (E-OD05), and Save writes no such entry again.
test('a Session raw LOSAT cache entry whose search context is a list is dropped', async () => {
  const session = syntheticConservationSession();
  const [, second] = session.losatCache.entries;
  second.searchContext = ['a'.repeat(64)];
  const result = await importSession(await sessionLoadEvent(session));
  assert.equal(result.status, 'ok');
  assert.deepEqual(result.repairs.map((repair) => [repair.code, repair.path]), [
    ['cache-entry-dropped', 'losatCache.entries[1]']
  ]);
  assert.equal(state.losatCache.value.has(second.key), false);
});

for (const version of [39, 40]) {
  test(`version-${version} import rejects duplicate raw LOSAT cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session([rawEntry('raw-key'), rawEntry('raw-key')], []),
        version
      ),
      sessionSchemaDiagnostic
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
      sessionSchemaDiagnostic
    );
  });

  for (const invalidKey of [undefined, '', null, 1]) {
    test(`version-${version} import rejects a nucleotide raw cache key of ${String(invalidKey)}`, () => {
      assert.throws(
        () => validateSessionLosatArtifacts(
          session([rawEntry(invalidKey)], []),
          version
        ),
        sessionSchemaDiagnostic
      );
    });
  }

  test(`version-${version} import rejects an empty protein raw cache key`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session([proteinRawEntry('')], []),
        version
      ),
      sessionSchemaDiagnostic
    );
  });

  for (const field of ['losatCache', 'losatDerivedCache']) {
    test(`version-${version} import rejects a non-object ${field}`, () => {
      const malformed = session([], []);
      malformed[field] = [];
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        sessionSchemaDiagnostic
      );
    });

    test(`version-${version} import rejects a non-array ${field}.entries`, () => {
      const malformed = session([], []);
      malformed[field] = { entries: null };
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        sessionSchemaDiagnostic
      );
    });
  }
  test(`version-${version} import rejects branch-internal raw schema 3`, () => {
    const malformed = session([{ ...proteinRawEntry('raw-key'), schema: 3 }], []);
    assert.throws(
      () => validateSessionLosatArtifacts(malformed, version),
      sessionSchemaDiagnostic
    );
  });

  test(`version-${version} import rejects branch-internal derived schema 2`, () => {
    const malformed = session([], [{ ...derivedEntry('derived-key'), schema: 2 }]);
    assert.throws(
      () => validateSessionLosatArtifacts(malformed, version),
      sessionSchemaDiagnostic
    );
  });
}

for (const unsupportedVersion of [34, 35, 36, 37, 38]) {
  test(`session import rejects branch-internal version ${unsupportedVersion}`, async () => {
    alerts.length = 0;
    const file = await sessionFileOf({
      format: 'gbdraw-session',
      version: unsupportedVersion
    });
    const event = { target: { files: [file], value: 'selected' } };
    const result = await importSession(event);
    assert.equal(result.status, 'error');
    assert.equal(result.error.code, 'INPUT_INVALID');
    // OV-1253: the refusal names the version and the versions this gbdraw reads.
    assert.deepEqual(result.error.context, {
      field: 'schema', reason: 'SESSION_VERSION_UNSUPPORTED', sessionVersion: unsupportedVersion,
      readSessionVersions: '27–33, 39–42, 44, and 46–51'
    });
    assert.match(result.error.summary, new RegExp(`Unsupported session version: ${unsupportedVersion}\\.`));
    assert.deepEqual(state.errorLog.value, result.error);
    assert.deepEqual(alerts, []);
    assert.equal(event.target.value, '');
  });
}

// OV-952: an unsafe own key anywhere in the file (here under `resources`) is
// refused with the Session-schema diagnostic, not as an unknown failure.
for (const key of ['__proto__', 'constructor', 'prototype']) {
  test(`session import refuses the unsafe key ${key} under resources with a classified diagnostic`, async () => {
    alerts.length = 0;
    const text = JSON.stringify(syntheticConservationSession())
      .replace('"resources":{', `"resources":{"${key}":{},`);
    assert.ok(text.includes(`"${key}":{}`));
    const event = await sessionLoadEventOfBytes(text);
    const result = await importSession(event);
    assert.equal(result.status, 'error');
    assert.deepEqual(
      [result.error.code, result.error.operation, result.error.stage, result.error.context],
      ['INPUT_INVALID', 'session-load', 'request-validation', { field: 'schema', reason: 'FIELDS' }]
    );
    assert.ok(!JSON.stringify(result.error).includes(key));
    assert.deepEqual(state.errorLog.value, result.error);
    assert.deepEqual(alerts, []);
  });
}

// OV-3179: a Session 31-39 Load reads the stored draft over its committed request's projection, and each projected
// ring takes the stored draft's genetic code of its own ring, paired by the one matcher, not by index. The document is
// the Web Session 39 of whitelist-tab-keyword.v39 with two TLOSATX rings in the shape the main cdd31013 Web writer
// gave them (the bytes of conservation-fasta.v39): each ring's table (`addText`: kind `conservation-blast-file`, the
// name `<genome>.circular_conservation.tlosatx.tsv` behind its resource ID, no original name), its genome in the
// request and in `webFiles.conservationLosatFastaSources` (`addFile`: kind `conservation-fasta-file`, the prefixed name
// and the original one), and a draft as that writer's Generate adoption left it: entry k at table k's place under
// table k's projected name (the projection drops the resource-ID prefix, then and now), with the source key of its
// genome from before the Generate.
const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url))));
const LOSAT_RINGS = [{ ring: 'ring-a', label: 'Ring A', code: 11 }, { ring: 'ring-b', label: 'Ring B', code: 4 }];
/** @param {(series: object[]) => object[]} editDraft What the user did to the draft's series after the Generate. */
const losatRingSession39 = (editDraft) => {
  const document = fixture('whitelist-tab-keyword.v39.gbdraw-session.json.gz');
  const { resources: ringBytes } = fixture('conservation-fasta.v39.gbdraw-session.json.gz');
  const table = ringBytes['conservation-blast-files-1'];
  const genome = ringBytes['conservation-fasta-files-1'];
  const originalNames = document.webFiles.resourceOriginalNames;
  const tableName = (ring) => `${ring}.circular_conservation.tlosatx.tsv`;
  const ids = LOSAT_RINGS.map(({ ring }, index) => {
    const tableId = `conservation-resolved-blast-${index + 1}`;
    document.resources[tableId] = { kind: 'conservation-blast-file', name: `${tableId}-${tableName(ring)}`,
      type: 'text/tab-separated-values', size: table.size, lastModified: 0, encoding: 'base64', data: table.data };
    const genomeIds = [`conservation-resolved-fasta-${index + 1}`, `conservation-losat-fasta-files-${index + 1}`];
    genomeIds.forEach((id) => {
      document.resources[id] = { ...genome, kind: 'conservation-fasta-file', name: `${id}-${ring}.fna` };
      originalNames[id] = `${ring}.fna`;
    });
    return { tableId, genomeId: genomeIds[0], losatId: genomeIds[1] };
  });
  Object.assign(document.renderRequest.diagramOptions, {
    conservationBlastFiles: ids.map(({ tableId }) => ({ resourceId: tableId, representation: 'file' })),
    conservationFastaFiles: ids.map(({ genomeId }) => ({ resourceId: genomeId, representation: 'file' })),
    conservationReference: 'subject',
    conservationLabels: LOSAT_RINGS.map(({ label }) => label),
    conservationColors: LOSAT_RINGS.map(() => '#D9EAF7'),
    conservationRingWidth: null,
    conservationRingGap: null
  });
  Object.assign(document.webFiles, {
    conservationBlastSource: 'losat-cache', conservationLosatFastaSources: ids.map(({ losatId }) => losatId)
  });
  const adopted = LOSAT_RINGS.map(({ ring, label, code }, index) => ({
    sourceKey: `${ring}.fna|${genome.size}|${genome.lastModified}|0`, fileName: tableName(ring),
    sourceIndex: index, label, color: '#d9eaf7', losat_gencode: code
  }));
  Object.assign(document.config.circularConservation, {
    enabled: true, source: 'losat', losat_program: 'tblastx', reference: 'subject',
    labels: LOSAT_RINGS.map(({ label }) => label).join(','), series: editDraft(adopted)
  });
  return document;
};
const loadedRingCodes = async (document) => {
  const result = await importSession(await sessionLoadEvent(document), {
    transformLegacyResultSvg, confirmLoadRepairs: async () => true
  });
  assert.equal(result.status, 'ok', JSON.stringify(result.error));
  return state.activeDrawing().circularConservation.series.map(({ label, losat_gencode: code }) => ({ label, code }));
};
const OWN_CODES = LOSAT_RINGS.map(({ label, code }) => ({ label, code }));

// Move series inward after the Generate: each entry keeps its table's name, so the matcher's name pass pairs it (its
// place would give the same pairing).
test('a Session 39 saved after a move of its LOSAT rings loads each ring with its own genetic code', async () => {
  assert.deepEqual(await loadedRingCodes(losatRingSession39((series) => [...series].reverse())), OWN_CODES,
    'each ring keeps its own genetic code');
});

// Move series inward, then Add comparison sequence: that writer cleared the derived tables and reconciled the series
// against the genome files, so each entry carries its genome's name, which no table has, and its genome's place,
// which is its table's place (the request lists the genomes in table order). The matcher's place pass decides; the
// added genome has no table and no ring.
test('a Session 39 saved after a move and Add comparison sequence loads each ring with its own genetic code', async () => {
  const document = losatRingSession39((series) => [
    ...[...series].reverse().map((entry) => ({ ...entry, fileName: entry.sourceKey.split('|')[0] })),
    { sourceKey: 'ring-c.fna|24|5|0', fileName: 'ring-c.fna', sourceIndex: 2, label: 'ring-c', color: '#d9eaf7', losat_gencode: 2 }
  ]);
  assert.deepEqual(await loadedRingCodes(document), OWN_CODES, 'each ring keeps its own genetic code');
});
