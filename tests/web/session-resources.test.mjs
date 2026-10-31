import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { File } from 'node:buffer';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { fixtureDrawing } from './helpers/drawing-state.mjs';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const repoRoot = process.cwd();
const sourceRoot = join(repoRoot, 'gbdraw', 'web', 'js');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-session-resources-'));
await cp(sourceRoot, join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const { buildSessionResources } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-resources.js'))
);

const makeFile = (text, name, options = {}) => new File([text], name, {
  type: options.type || 'text/plain',
  lastModified: options.lastModified ?? 1
});
const base64 = (text) => Buffer.from(text).toString('base64');
const activeText = 'byte-identical';
const activeFile = makeFile(activeText, 'circular-display.gb', {
  type: 'text/x-genbank',
  lastModified: 101
});
let activeReads = 0;
const originalArrayBuffer = activeFile.arrayBuffer.bind(activeFile);
activeFile.arrayBuffer = async () => {
  activeReads += 1;
  return originalArrayBuffer();
};
const separatelyConstructedDuplicate = makeFile(
  activeText,
  'linear-display.gb',
  { type: 'application/genbank', lastModified: 202 }
);
const sameMetadataDifferentA = makeFile('different-A', 'same.tsv', {
  lastModified: 303
});
const sameMetadataDifferentB = makeFile('different-B', 'same.tsv', {
  lastModified: 303
});

const committed = {
  renderRequest: {
    schema: CANONICAL_REQUEST_SCHEMA,
    mode: 'circular',
    records: [
      {
        recordKey: 'record-1',
        source: { kind: 'genbank', resourceId: 'record-source' },
        presentation: { label: 'record-source' }
      },
      {
        recordKey: 'record-2',
        source: {
          kind: 'gffFasta',
          gffResourceId: 'record-gff',
          fastaResourceId: 'record-fasta'
        },
        presentation: { label: 'GFF record' }
      }
    ],
    diagramOptions: {
      trackSlots: [{
        id: 'record-source',
        renderer: 'annotations',
        params: {
          anchor_slot: 'record-source',
          label: 'record-source'
        }
      }],
      annotations: {
        sets: [{
          id: 'record-source',
          annotations: [{ id: 'record-source', label: 'record-source' }]
        }]
      }
    },
    layout: {},
    comparisons: [],
    output: { prefix: 'record-source', formats: ['svg'], overwrite: false }
  },
  resources: {
    'record-source': {
      kind: 'genbank',
      name: 'record-source-active.gb',
      type: 'text/x-genbank',
      size: Buffer.byteLength(activeText),
      lastModified: 0,
      encoding: 'base64',
      data: base64(activeText)
    },
    'record-gff': {
      kind: 'gff3',
      name: 'record-gff.gff3',
      type: 'text/plain',
      size: Buffer.byteLength('##gff-version 3\n'),
      lastModified: 0,
      encoding: 'base64',
      data: base64('##gff-version 3\n')
    },
    'record-fasta': {
      kind: 'fasta',
      name: 'record-fasta.fa',
      type: 'text/plain',
      size: Buffer.byteLength('>record\nATGC\n'),
      lastModified: 0,
      encoding: 'base64',
      data: base64('>record\nATGC\n')
    },
    'removed-inactive-resource': {
      kind: 'web-file',
      name: 'removed-inactive.txt',
      type: 'text/plain',
      size: 24,
      lastModified: 0,
      encoding: 'base64',
      data: base64('removed-inactive-content')
    }
  },
  webFiles: {
    resourceOriginalNames: {
      'record-source': 'original-active.gb',
      'removed-inactive-resource': 'removed-inactive.txt'
    },
    conservationLosatFastaSources: ['removed-inactive-resource']
  }
};

const state = {
  files: {
    c_gb: activeFile,
    c_gff: sameMetadataDifferentA,
    c_fasta: sameMetadataDifferentB,
    c_depth: [[makeFile('depth-1', 'depth-1.tsv')]],
    c_conservation_blasts: [makeFile('blast', 'blast.tsv')],
    c_conservation_blasts_source: 'losat-cache',
    c_conservation_fastas: [makeFile('>subject\nATGC\n', 'subject.fa')],
    c_conservation_sequence_sources: [makeFile('>source\nATGC\n', 'source.fa')],
    d_color: activeFile,
    t_color: makeFile('CDS\t#ffffff\n', 'specific.tsv'),
    blacklist: makeFile('hypothetical', 'blacklist.txt'),
    whitelist: makeFile('CDS\tgene\tfoo\n', 'whitelist.tsv'),
    qualifier_priority: makeFile('CDS\tgene\n', 'priority.tsv'),
    linearCanonicalComparisons: [{
      kind: 'collinearityResult',
      encoding: 'canonicalJson',
      valueKind: 'result',
      file: makeFile('{"blocks":[]}', 'collinearity.json')
    }]
  },
  linearSeqs: [{
    uid: 'linear-uid-1',
    gb: separatelyConstructedDuplicate,
    gff: null,
    fasta: null,
    depth: [makeFile('linear-depth', 'linear-depth.tsv')],
    losat_gencode: 11,
    definition: 'Linear record',
    record_subtitle: 'Subtitle',
    file_definition: 'Default organism',
    file_subtitle: 'Default subtitle',
    inferred_definition: '<i>Escherichia coli</i> K-12',
    region_record_id: '#1',
    region_start: 2,
    region_end: 9,
    region_reverse: true
  }],
  linearComparisonPlan: {
    mode: 'selected',
    defaultSource: 'losat',
    edges: [{
      id: 'comparison-uid-1',
      queryUid: 'linear-uid-1',
      subjectUid: 'linear-uid-2',
      included: true,
      fileActive: true,
      losatFilenameActive: false,
      source: 'upload',
      file: makeFile('uploaded-comparison', 'comparison.tsv'),
      losatFilename: ''
    }]
  }
};

const built = await buildSessionResources(state, committed, fixtureDrawing(state));
const bindings = built.webFiles.bindings;
assert.equal(bindings.schema, 2);
assert.equal(activeReads, 1, 'one File object used in several roles is read once');
assert.equal(
  bindings.c_gb.resourceId,
  bindings.d_color.resourceId,
  'one File object used in several roles shares payload bytes'
);
assert.equal(
  bindings.c_gb.resourceId,
  bindings.linearSeqs[0].gb.resourceId,
  'separate byte-identical File objects share payload bytes'
);
assert.notEqual(
  bindings.c_gff.resourceId,
  bindings.c_fasta.resourceId,
  'same metadata with different bytes must not deduplicate'
);
assert.deepEqual(bindings.c_gb, {
  resourceId: bindings.c_gb.resourceId,
  name: 'circular-display.gb',
  type: 'text/x-genbank',
  lastModified: 101
});
assert.deepEqual(bindings.linearSeqs[0].gb, {
  resourceId: bindings.c_gb.resourceId,
  name: 'linear-display.gb',
  type: 'application/genbank',
  lastModified: 202
});
assert.equal(bindings.linearSeqs[0].uid, 'linear-uid-1');
assert.equal(bindings.linearSeqs[0].file_definition, 'Default organism');
assert.equal(bindings.linearSeqs[0].file_subtitle, 'Default subtitle');
// IN-06 (D-12): the current writer keeps each record's inferred definition.
assert.equal(bindings.linearSeqs[0].inferred_definition, '<i>Escherichia coli</i> K-12');
assert.equal(bindings.linearComparisons[0].id, 'comparison-uid-1');
assert.deepEqual(Object.keys(bindings.linearComparisons[0]), ['id', 'file']);
assert.equal(Object.hasOwn(bindings, 'linearCanonicalComparisons'), false);
assert.equal(
  built.renderRequest.records[0].source.resourceId,
  bindings.c_gb.resourceId
);
assert.equal(
  built.renderRequest.records[1].source.gffResourceId in built.resources,
  true,
  'GFF resource references are retained and remapped'
);
assert.equal(
  built.renderRequest.records[1].source.fastaResourceId in built.resources,
  true,
  'FASTA resource references are retained and remapped'
);
assert.notEqual(
  built.renderRequest.records[1].source.gffResourceId,
  'record-gff'
);
assert.notEqual(
  built.renderRequest.records[1].source.fastaResourceId,
  'record-fasta'
);
assert.equal(
  built.renderRequest.records[0].presentation.label,
  'record-source',
  'record labels that happen to match a prior resource ID are not rewritten'
);
assert.equal(
  built.renderRequest.diagramOptions.trackSlots[0].id,
  'record-source',
  'slot IDs that happen to match a prior resource ID are not rewritten'
);
assert.equal(
  built.renderRequest.diagramOptions.trackSlots[0].params.anchor_slot,
  'record-source',
  'annotation anchors that happen to match a prior resource ID are not rewritten'
);
assert.equal(
  built.renderRequest.diagramOptions.annotations.sets[0].annotations[0].label,
  'record-source',
  'annotation labels that happen to match a prior resource ID are not rewritten'
);
assert.equal(
  built.renderRequest.output.prefix,
  'record-source',
  'output prefixes that happen to match a prior resource ID are not rewritten'
);
assert.deepEqual(
  Object.keys(built.webFiles.resourceOriginalNames),
  [bindings.c_gb.resourceId]
);
assert.deepEqual(built.webFiles.conservationLosatFastaSources, [null]);
assert.equal(
  Object.values(built.resources).some(
    (resource) => Buffer.from(resource.data, 'base64').toString('utf8')
      === 'removed-inactive-content'
  ),
  false,
  'unreferenced resources from a prior session are not retained after the file is removed'
);
assert.ok(
  Object.keys(built.resources).every((resourceId) => /^resource-\d{4}$/.test(resourceId)),
  'resource IDs stay opaque and do not expose digests or file names'
);
assert.ok(
  Object.entries(built.resources).every(
    ([resourceId, resource]) => resource.name.startsWith(`${resourceId}-`)
  )
);

const dormantComparisonFile = makeFile(
  'dormant-uploaded-comparison',
  'dormant-comparison.tsv'
);
let dormantComparisonReads = 0;
const dormantComparisonArrayBuffer = dormantComparisonFile.arrayBuffer.bind(
  dormantComparisonFile
);
dormantComparisonFile.arrayBuffer = async () => {
  dormantComparisonReads += 1;
  return dormantComparisonArrayBuffer();
};
const dormantSaveState = {
  ...state,
  linearComparisonPlan: {
    mode: 'none',
    defaultSource: 'losat',
    edges: [{
      id: 'dormant-comparison-uid',
      queryUid: 'linear-uid-1',
      subjectUid: 'linear-uid-2',
      included: false,
      fileActive: false,
      losatFilenameActive: false,
      source: 'upload',
      file: dormantComparisonFile,
      losatFilename: 'retained.raw.tsv'
    }]
  }
};
const dormantCommitted = structuredClone(committed);
dormantCommitted.renderRequest.mode = 'linear';
dormantCommitted.renderRequest.comparisons = [];
const dormantSaved = await buildSessionResources(
  dormantSaveState,
  dormantCommitted,
  fixtureDrawing(dormantSaveState)
);
assert.equal(
  dormantComparisonReads,
  1,
  'Save Session reads an inactive comparison file so its editable binding survives'
);
assert.deepEqual(dormantSaved.renderRequest.comparisons, []);
assert.equal(
  dormantSaved.webFiles.bindings.linearComparisons[0].id,
  'dormant-comparison-uid'
);
const dormantResourceId =
  dormantSaved.webFiles.bindings.linearComparisons[0].file.resourceId;
assert.equal(
  Buffer.from(dormantSaved.resources[dormantResourceId].data, 'base64').toString('utf8'),
  'dormant-uploaded-comparison'
);

// C1 I10/I11/I13/I22/I23: payload identity is independent of draft metadata.
const backing = await import(pathToFileURL(join(tempRoot, 'js/services/session-resource-backing.js')));
const authority = await import(pathToFileURL(join(tempRoot, 'js/services/session-authority.js')));
const { readFile } = await import('node:fs/promises');
const { createHash } = await import('node:crypto');
const frozenBytes = await readFile(join(repoRoot, 'tests/fixtures/sessions/single.v41-bindings1.json'));
const frozenProvenance = JSON.parse(await readFile(join(repoRoot, 'tests/fixtures/sessions/single.v41-bindings1.provenance.json')));
assert.equal(createHash('sha256').update(frozenBytes).digest('hex'), frozenProvenance.sha256);
const frozen = JSON.parse(frozenBytes);
const frozenLeaf = frozen.webFiles.bindings.c_gb;
const frozenTable = backing.adoptCurrentSessionResources(frozen.resources);
const oldFile = backing.createSessionResourceFileView(frozenTable, frozenLeaf.resourceId, frozenLeaf);
const promoted = await buildSessionResources({ files: { c_gb: oldFile } }, authority.adoptRuntimeCanonicalSession(frozen));
assert.equal(promoted.webFiles.bindings.schema, 2);
assert.deepEqual(promoted.webFiles.bindings.c_gb, frozenLeaf);

const descriptor = (text, name = 'same.gb') => ({ kind: 'genbank', name, type: 'text/plain',
  encoding: 'base64', size: Buffer.byteLength(text), data: base64(text), lastModified: 0 });
for (const adopted of [false, true]) {
  const committed = { renderRequest: { schema: CANONICAL_REQUEST_SCHEMA, records: [{ source: { resourceId: 'occupied' } }] },
    resources: { occupied: descriptor('committed\n', 'committed.gb'), same: descriptor('A\n') } };
  const parts = { occupied: descriptor('A\n'), unused: descriptor('B\n', 'unique.gb') };
  const table = backing.adoptCurrentSessionResources(parts);
  const components = [
    { resourceId: 'unused', name: '', type: '', lastModified: 0.5 },
    { resourceId: 'occupied', name: 'one.gb', type: 'text/plain', lastModified: 1 },
    { resourceId: 'occupied', name: 'repeat.gb', type: 'different', lastModified: 2 }
  ];
  const file = backing.createCombinedSessionResourceFileView(table, components, {
    name: '', type: '', lastModified: 0.25
  });
  const state = { files: { c_gb: file, c_fasta: makeFile('A\n', 'independent.fa') } };
  const canonical = adopted ? authority.adoptRuntimeCanonicalSession(committed) : committed;
  const first = await buildSessionResources(state, canonical, fixtureDrawing(state));
  const binding = first.webFiles.bindings.c_gb;
  assert.equal(binding.kind, 'composite');
  assert.equal(binding.name, '');
  assert.equal(binding.components[0].name, '');
  assert.deepEqual(binding.components.map(c => first.resources[c.resourceId].data),
    ['B\n', 'A\n', 'A\n'].map(base64));
  assert.equal(binding.components[1].resourceId, binding.components[2].resourceId);
  assert.equal(binding.components[1].resourceId, first.webFiles.bindings.c_fasta.resourceId);
  if (adopted) {
    assert.equal(binding.components[1].resourceId, 'same');
    assert.equal(binding.components[0].resourceId, 'unused');
    assert.equal(first.renderRequest, committed.renderRequest);
  }
  const loaded = backing.createCombinedSessionResourceFileView(
    backing.adoptCurrentSessionResources(first.resources), binding.components, binding
  );
  const second = await buildSessionResources({ files: { c_gb: loaded } }, authority.adoptRuntimeCanonicalSession(first));
  assert.deepEqual(second.webFiles.bindings.c_gb, binding);
  state.files.c_gb = makeFile('native replacement\n', 'replacement.gb');
  const replaced = await buildSessionResources(state, canonical, fixtureDrawing(state));
  assert.equal(replaced.webFiles.bindings.c_gb.kind, undefined);
  assert.equal(replaced.webFiles.bindings.c_gb.components, undefined);
  assert.equal(replaced.resources[replaced.webFiles.bindings.c_gb.resourceId].data, base64('native replacement\n'));
}

// E1: the other mode's committed request (`otherModeResult`) shares the
// Session's one resource table. Equal bytes in both requests share one
// resource; a positional ID that names different bytes in the two requests
// gets a second resource, and each request's references are rewritten.
{
  const { adoptRuntimeCanonicalSession } = await import(
    pathToFileURL(join(tempRoot, 'js', 'services', 'session-authority.js'))
  );
  const genbank = (text, name) => ({
    kind: 'genbank', name, type: 'text/x-genbank', size: Buffer.byteLength(text),
    lastModified: 0, encoding: 'base64', data: base64(text)
  });
  const request = (mode) => ({
    schema: 9, mode, diagramOptions: {}, layout: {}, comparisons: [],
    records: ['record-1-genbank', 'record-2-genbank'].map((resourceId, index) => ({
      recordKey: `record-${index + 1}`, source: { kind: 'genbank', resourceId }
    })),
    output: { prefix: mode, formats: ['svg'], overwrite: false }
  });
  const circular = {
    renderRequest: request('circular'),
    resources: { 'record-1-genbank': genbank('CIRCULAR', 'record-1.gb'), 'record-2-genbank': genbank('SHARED', 'record-2.gb') },
    webFiles: { resourceOriginalNames: { 'record-1-genbank': 'circular.gb' } }
  };
  const linear = {
    renderRequest: request('linear'),
    resources: { 'record-1-genbank': genbank('LINEAR', 'record-1.gb'), 'record-2-genbank': genbank('SHARED', 'record-2.gb') },
    webFiles: { resourceOriginalNames: { 'record-1-genbank': 'linear.gb' }, linearRecordMetadata: [{ recordKey: 'record-1' }] }
  };
  const noInputs = { files: {}, linearSeqs: [], linearComparisonPlan: { edges: [] } };
  const recordResource = (built, renderRequest, index) => built.resources[renderRequest.records[index].source.resourceId];
  const text = (descriptor) => Buffer.from(descriptor.data, 'base64').toString();

  const both = await buildSessionResources(noInputs, circular, fixtureDrawing(noInputs), linear);
  assert.equal(text(recordResource(both, both.renderRequest, 0)), 'CIRCULAR');
  assert.equal(text(recordResource(both, both.otherRenderRequest, 0)), 'LINEAR');
  assert.notEqual(both.renderRequest.records[0].source.resourceId, both.otherRenderRequest.records[0].source.resourceId);
  assert.equal(both.renderRequest.records[1].source.resourceId, both.otherRenderRequest.records[1].source.resourceId,
    'equal bytes in both modes share one resource');
  assert.equal(Object.keys(both.resources).length, 3);
  assert.equal(both.otherRenderRequest.mode, 'linear');
  assert.deepEqual(both.webFiles.linearRecordMetadata, [{ recordKey: 'record-1' }]);
  assert.equal(both.webFiles.resourceOriginalNames[both.renderRequest.records[0].source.resourceId], 'circular.gb');
  assert.equal(both.webFiles.resourceOriginalNames[both.otherRenderRequest.records[0].source.resourceId], 'linear.gb');

  // Both requests adopted from one loaded Session's table keep their IDs.
  const table = {
    ...circular.resources,
    'linear-record-1-genbank': genbank('LINEAR', 'linear-record-1.gb')
  };
  const linearRequest = request('linear');
  linearRequest.records[0].source.resourceId = 'linear-record-1-genbank';
  const adopted = await buildSessionResources(
    noInputs,
    adoptRuntimeCanonicalSession({ renderRequest: circular.renderRequest, resources: table, webFiles: {} }),
    fixtureDrawing(noInputs),
    adoptRuntimeCanonicalSession({ renderRequest: linearRequest, resources: table, webFiles: {} })
  );
  assert.strictEqual(adopted.renderRequest, circular.renderRequest);
  assert.strictEqual(adopted.otherRenderRequest, linearRequest);
  assert.deepEqual(Object.keys(adopted.resources).sort(), Object.keys(table).sort());
  assert.strictEqual(adopted.resources['linear-record-1-genbank'], table['linear-record-1-genbank']);
}

// E1 (review m4): a Session loaded with a Result of each mode, whose Linear
// Result is then regenerated from a replaced file, writes the new Linear bytes
// and not the replaced ones: only resources the requests and bindings name.
{
  const { adoptRuntimeCanonicalSession } = await import(
    pathToFileURL(join(tempRoot, 'js', 'services', 'session-authority.js'))
  );
  const genbank = (text, name) => ({
    kind: 'genbank', name, type: 'text/x-genbank', size: Buffer.byteLength(text),
    lastModified: 0, encoding: 'base64', data: base64(text)
  });
  const request = (mode, ids) => ({
    schema: 9, mode, diagramOptions: {}, layout: {}, comparisons: [],
    records: ids.map((resourceId, index) => ({ recordKey: `record-${index + 1}`, source: { kind: 'genbank', resourceId } })),
    output: { prefix: mode, formats: ['svg'], overwrite: false }
  });
  const table = {
    'record-1-genbank': genbank('CIRCULAR', 'record-1.gb'),
    'linear-record-1-genbank': genbank('OLD LINEAR', 'linear-record-1.gb')
  };
  const loadedCircular = adoptRuntimeCanonicalSession({
    renderRequest: request('circular', ['record-1-genbank']), resources: table, webFiles: {}
  });
  const freshLinear = {
    renderRequest: request('linear', ['record-1-genbank']),
    resources: { 'record-1-genbank': genbank('NEW LINEAR', 'new.gb') },
    webFiles: {}
  };
  const draft = { files: {}, linearSeqs: [], linearComparisonPlan: { edges: [] } };
  const built = await buildSessionResources(draft, loadedCircular, fixtureDrawing(draft), freshLinear);
  const written = Object.values(built.resources).map((descriptor) => Buffer.from(descriptor.data, 'base64').toString()).sort();
  assert.deepEqual(written, ['CIRCULAR', 'NEW LINEAR']);
}

// E1 (REVIEW-2 W): a Circular Result kept in `otherModeResult` while Linear is
// shown keeps the request metadata only a Circular commit writes: the input's
// original file name and the LOSAT-cache conservation marker. They survive a
// Save, the Load that adopts the two committed Sessions, and the next Save.
{
  const { adoptRuntimeCanonicalSession } = await import(
    pathToFileURL(join(tempRoot, 'js', 'services', 'session-authority.js'))
  );
  const genbank = (text, name) => ({
    kind: 'genbank', name, type: 'text/x-genbank', size: Buffer.byteLength(text),
    lastModified: 0, encoding: 'base64', data: base64(text)
  });
  const request = (mode, id) => ({
    schema: 9, mode, diagramOptions: {}, layout: {}, comparisons: [],
    records: [{ recordKey: 'record-1', source: { kind: 'genbank', resourceId: id } }],
    output: { prefix: mode, formats: ['svg'], overwrite: false }
  });
  const shownLinear = {
    renderRequest: request('linear', 'record-1-genbank'),
    resources: { 'record-1-genbank': genbank('LINEAR', 'linear.gb') },
    webFiles: { linearRecordMetadata: [{ recordKey: 'record-1', definition: 'Linear record' }] }
  };
  const keptCircular = {
    renderRequest: request('circular', 'record-1-genbank'),
    resources: { 'record-1-genbank': genbank('CIRCULAR', 'circular.gb') },
    webFiles: { circularInputOriginalName: 'genome.gbk', conservationBlastSource: 'losat-cache' }
  };
  const draft = { files: {}, linearSeqs: [], linearComparisonPlan: { edges: [] } };
  const saved = await buildSessionResources(draft, shownLinear, fixtureDrawing(draft), keptCircular);
  assert.equal(saved.webFiles.circularInputOriginalName, 'genome.gbk');
  assert.equal(saved.webFiles.conservationBlastSource, 'losat-cache');
  assert.deepEqual(saved.webFiles.linearRecordMetadata, shownLinear.webFiles.linearRecordMetadata);
  const loaded = (renderRequest) => adoptRuntimeCanonicalSession({
    renderRequest, resources: saved.resources, webFiles: saved.webFiles
  });
  const resaved = await buildSessionResources(draft, loaded(saved.renderRequest), fixtureDrawing(draft), loaded(saved.otherRenderRequest));
  assert.equal(resaved.webFiles.circularInputOriginalName, 'genome.gbk');
  assert.equal(resaved.webFiles.conservationBlastSource, 'losat-cache');
}

// D-220, D-228 (Owner 2026-10-30): a committed request that names two
// resources with the same bytes, kind and original name (one Depth file added
// to two series) is saved with one resource; every reference names it, and the
// Load of the saved bindings gives each series the bytes. Equal names with
// different bytes stay two.
{
  const depth = (text, name) => ({
    kind: 'depth-track-file', name, type: 'text/tab-separated-values', size: Buffer.byteLength(text),
    lastModified: 0, encoding: 'base64', data: base64(text)
  });
  const track = (index, ids) => ({ id: `depth-${index}`, source: ids.map((resourceId) => ({ resourceId })) });
  const renderRequest = {
    schema: CANONICAL_REQUEST_SCHEMA, mode: 'linear', layout: {}, comparisons: [],
    records: [1, 2].map((n) => ({ recordKey: `record-${n}`, source: { kind: 'genbank', resourceId: `record-${n}-genbank` } })),
    diagramOptions: { depthTracks: [
      track(1, ['depth-tracks-1-source-record-1', 'depth-tracks-1-source-record-2']),
      track(2, ['depth-tracks-2-source-record-1', 'depth-tracks-2-source-record-2'])
    ] },
    output: { prefix: 'x', formats: ['svg'], overwrite: false }
  };
  const resources = {
    'record-1-genbank': { ...depth('LOCUS one\n', 'record-1-genbank-one.gb'), kind: 'genbank' },
    'record-2-genbank': { ...depth('LOCUS two\n', 'record-2-genbank-two.gb'), kind: 'genbank' },
    'depth-tracks-1-source-record-1': depth('chr\t1\t5\n', 'depth-tracks-1-source-record-1-file1.tsv'),
    'depth-tracks-2-source-record-1': depth('chr\t1\t5\n', 'depth-tracks-2-source-record-1-file1.tsv'),
    'depth-tracks-1-source-record-2': depth('chr\t1\t6\n', 'depth-tracks-1-source-record-2-same.tsv'),
    'depth-tracks-2-source-record-2': depth('chr\t1\t7\n', 'depth-tracks-2-source-record-2-same.tsv')
  };
  const webFiles = { resourceOriginalNames: {
    'depth-tracks-1-source-record-1': 'file1.tsv', 'depth-tracks-2-source-record-1': 'file1.tsv',
    'depth-tracks-1-source-record-2': 'same.tsv', 'depth-tracks-2-source-record-2': 'same.tsv'
  } };
  const table = backing.adoptCurrentSessionResources(resources);
  const view = (id) => backing.createSessionResourceFileView(table, id);
  const draft = {
    files: {},
    linearSeqs: [1, 2].map((n) => ({ uid: `uid-${n}`, gb: view(`record-${n}-genbank`),
      depth: [view(`depth-tracks-1-source-record-${n}`), view(`depth-tracks-2-source-record-${n}`)] })),
    linearComparisonPlan: { edges: [] }
  };
  const committedRequest = authority.adoptRuntimeCanonicalSession({ renderRequest, resources, webFiles });
  const saved = await buildSessionResources(draft, committedRequest, fixtureDrawing(draft));
  const tracks = saved.renderRequest.diagramOptions.depthTracks;
  const kept = 'depth-tracks-1-source-record-1';
  assert.deepEqual(tracks.map((t) => t.source[0].resourceId), [kept, kept], 'one resource for equal bytes');
  assert.equal(Object.hasOwn(saved.resources, 'depth-tracks-2-source-record-1'), false);
  assert.strictEqual(saved.resources[kept], resources[kept]);
  assert.deepEqual(saved.webFiles.bindings.linearSeqs[0].depth.map((b) => b.resourceId), [kept, kept]);
  assert.equal(saved.webFiles.resourceOriginalNames[kept], 'file1.tsv', 'the kept resource keeps its own name');
  assert.deepEqual(tracks.map((t) => t.source[1].resourceId),
    ['depth-tracks-1-source-record-2', 'depth-tracks-2-source-record-2'], 'equal names, different bytes: two');
  assert.deepEqual(Object.keys(saved.resources).sort(), Object.keys(resources).filter((id) => id !== 'depth-tracks-2-source-record-1').sort());
  assert.notStrictEqual(saved.renderRequest, renderRequest);
  assert.equal(renderRequest.diagramOptions.depthTracks[1].source[0].resourceId, 'depth-tracks-2-source-record-1', 'the committed request is not edited');
  // Load: each saved series binding reads its own bytes back.
  const loadedTable = backing.adoptCurrentSessionResources(saved.resources);
  const loadedDepth = await Promise.all(saved.webFiles.bindings.linearSeqs.flatMap((seq) => seq.depth).map(async (binding) => (
    Buffer.from(await backing.readSessionResourceBytes(backing.createSessionResourceFileView(loadedTable, binding.resourceId, binding))).toString()
  )));
  assert.deepEqual(loadedDepth, ['chr\t1\t5\n', 'chr\t1\t5\n', 'chr\t1\t6\n', 'chr\t1\t7\n']);
  // A second Save of the loaded Session writes the same table and request.
  const reloaded = backing.adoptCurrentSessionResources(saved.resources);
  const redraft = { ...draft, linearSeqs: saved.webFiles.bindings.linearSeqs.map((seq) => ({ uid: seq.uid,
    gb: backing.createSessionResourceFileView(reloaded, seq.gb.resourceId, seq.gb),
    depth: seq.depth.map((b) => backing.createSessionResourceFileView(reloaded, b.resourceId, b)) })) };
  const resaved = await buildSessionResources(redraft, authority.adoptRuntimeCanonicalSession({
    renderRequest: saved.renderRequest, resources: saved.resources, webFiles: saved.webFiles
  }), fixtureDrawing(redraft));
  assert.deepEqual(Object.keys(resaved.resources).sort(), Object.keys(saved.resources).sort());
  assert.strictEqual(resaved.renderRequest, saved.renderRequest);
}

// D-228: equal bytes merge only with an equal kind and an equal original name
// (a missing name equals only a missing name), looked up per kind and name, so
// of X, Y, Y with equal bytes the two Y copies merge.
{
  const descriptorOf = (kind, id) => ({ kind, name: `${id}.dat`, type: 'text/plain', size: 5,
    lastModified: 0, encoding: 'base64', data: base64('SAME\n') });
  const save = async (entries, names) => {
    const ids = entries.map(([id]) => id);
    const renderRequest = { schema: CANONICAL_REQUEST_SCHEMA, mode: 'circular', layout: {}, comparisons: [],
      records: [], diagramOptions: { depthTracks: ids.map((resourceId, index) => ({ id: `t${index}`, source: { resourceId } })) },
      output: { prefix: 'x', formats: ['svg'], overwrite: false } };
    const resources = Object.fromEntries(entries.map(([id, kind]) => [id, descriptorOf(kind, id)]));
    const draft = { files: {}, linearSeqs: [], linearComparisonPlan: { edges: [] } };
    const saved = await buildSessionResources(draft, authority.adoptRuntimeCanonicalSession({
      renderRequest, resources, webFiles: { resourceOriginalNames: names }
    }), fixtureDrawing(draft));
    return saved.renderRequest.diagramOptions.depthTracks.map((track) => track.source.resourceId);
  };
  const depthKind = 'depth-track-file';
  assert.deepEqual(await save([['a', depthKind], ['b', depthKind]], { a: 'a.tsv', b: 'b.tsv' }), ['a', 'b'],
    'different names stay two');
  assert.deepEqual(await save([['a', depthKind], ['b', depthKind]], { a: 'a.tsv' }), ['a', 'b'],
    'a missing name differs from a name');
  assert.deepEqual(await save([['a', depthKind], ['b', depthKind]], {}), ['a', 'a'],
    'two missing names are equal');
  assert.deepEqual(await save([['x', 'fasta'], ['y1', depthKind], ['y2', depthKind]], { x: 's.tsv', y1: 's.tsv', y2: 's.tsv' }),
    ['x', 'y1', 'y1'], 'X, Y, Y: the two Y copies merge');
}
