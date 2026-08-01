import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { File } from 'node:buffer';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

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
    schema: 5,
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

const built = await buildSessionResources(state, committed);
const bindings = built.webFiles.bindings;
assert.equal(bindings.schema, 1);
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
