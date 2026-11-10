import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { createHash, webcrypto } from 'node:crypto';
import test from 'node:test';
import { bindingDrawings, boundState } from './helpers/input-files.mjs';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric(metric) {
    metrics.push({ ...metric });
  }
};

const {
  adoptCurrentSessionResources,
  createSessionResourceFileView,
  readSessionResourceBytes,
  readSessionResourceText,
  validateWebFileBindings
} = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const {
  holdReadBytes,
  readFileBytes,
  readFileContentFacts,
  readFileText
} = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { adoptRuntimeCanonicalSession } = await import(
  '../../gbdraw/web/js/services/session-authority.js'
);
const { createInputFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const { buildSessionResources: buildBoundSessionResources } = await import(
  '../../gbdraw/web/js/services/session-resources.js'
);
const { writtenResourcesWith } = await import('./helpers/written-resources.mjs');
// What Save writes for each resource, read as base64 descriptors.
const encodeSessionResources = writtenResourcesWith(await import('../../gbdraw/web/js/services/session-resource-backing.js'));
// A fixture state registers its input files in a store and names them by ref, as the app does;
// its drawings are those of the bound state, so both name the same refs. Save then encodes the
// table in the Session file writer.
const buildSessionResources = async (state, committed, drawings = bindingDrawings(state), ...rest) => {
  const bound = state.inputFiles ? state : boundState(state);
  const built = await buildBoundSessionResources(bound, committed,
    state.inputFiles ? drawings : bindingDrawings(bound), ...rest);
  return { ...built, resources: await encodeSessionResources(built.resources) };
};

const encodedDescriptor = (name, text, overrides = {}) => {
  const bytes = new TextEncoder().encode(text);
  return {
    kind: 'genbank',
    name,
    type: 'text/plain',
    size: bytes.byteLength,
    lastModified: 7,
    encoding: 'base64',
    data: Buffer.from(bytes).toString('base64'),
    ...overrides
  };
};

const metricTotal = (name, resourceId = null) => metrics
  .filter((metric) => (
    metric.name === name
    && (resourceId === null || metric.resourceId === resourceId)
  ))
  .reduce((total, metric) => total + metric.value, 0);

const resetMetrics = () => metrics.splice(0);

test('current-session file views materialize one requested resource once', async () => {
  const alphaText = 'alpha resource\n';
  const checksum = createHash('sha256').update(alphaText).digest('hex');
  const resources = {
    alpha: encodedDescriptor('alpha.gbk', alphaText, {
      checksum: `sha256:${checksum}`
    }),
    beta: encodedDescriptor('beta.gbk', 'beta resource\n')
  };
  resetMetrics();
  const table = adoptCurrentSessionResources(resources);
  const alpha = createSessionResourceFileView(table, 'alpha', {
    name: 'display-alpha.gbk',
    type: 'application/genbank',
    lastModified: 99
  });
  const beta = createSessionResourceFileView(table, 'beta');

  assert.deepEqual(alpha, {
    name: 'display-alpha.gbk',
    type: 'application/genbank',
    size: new TextEncoder().encode(alphaText).byteLength,
    lastModified: 99
  });
  assert.equal(Object.isFrozen(alpha), true);
  assert.equal(Object.hasOwn(alpha, 'data'), false);
  assert.equal(Object.hasOwn(alpha, 'resourceId'), false);
  assert.equal(metrics.length, 0, 'adoption and metadata projection stay lazy');

  const directBytes = readSessionResourceBytes(alpha);
  assert.strictEqual(
    readSessionResourceBytes(alpha),
    directBytes,
    'concurrent backing reads share one Promise'
  );
  const genericBytes = readFileBytes(alpha);
  assert.strictEqual(
    genericBytes,
    directBytes,
    'the generic cache retains the backing Promise rather than a copied byte array'
  );
  assert.strictEqual(
    readFileBytes(alpha),
    genericBytes,
    'concurrent generic reads share the cached Promise'
  );
  const [directValue, bytesA, bytesB, textA, textB, directText] = await Promise.all([
    directBytes,
    genericBytes,
    readFileBytes(alpha),
    readFileText(alpha),
    readFileText(alpha),
    readSessionResourceText(alpha)
  ]);

  assert.strictEqual(bytesA, directValue);
  assert.strictEqual(bytesA, bytesB);
  assert.strictEqual(
    await readFileBytes(alpha),
    directValue,
    'the generic cache does not retain a second Uint8Array'
  );
  assert.equal(new TextDecoder().decode(bytesA), alphaText);
  assert.equal(textA, alphaText);
  assert.equal(textB, alphaText);
  assert.equal(directText, alphaText);
  assert.equal(metricTotal('base64DecodeCount', 'alpha'), 1);
  assert.equal(metricTotal('decodedByteCount', 'alpha'), alpha.size);
  assert.equal(metricTotal('resourceByteReadCount', 'alpha'), 1);
  assert.equal(metricTotal('resourceTextReadCount', 'alpha'), 1);
  assert.equal(metricTotal('base64DecodeCount', 'beta'), 0);
  assert.equal(metricTotal('resourceByteReadCount', 'beta'), 0);
  assert.equal(await readFileText(beta), 'beta resource\n');
  assert.equal(metricTotal('base64DecodeCount', 'beta'), 1);
});

test('lazy resource failures identify the resource and remain cached', async () => {
  const cases = [
    {
      id: 'invalid-base64',
      descriptor: encodedDescriptor('invalid.gbk', 'unused', { data: '%%%' }),
      message: /invalid-base64 \(invalid\.gbk\) contains invalid encoded data/
    },
    {
      id: 'wrong-size',
      descriptor: encodedDescriptor('wrong-size.gbk', 'size', { size: 99 }),
      message: /wrong-size \(wrong-size\.gbk\) byte size does not match/
    },
    {
      id: 'wrong-checksum',
      descriptor: encodedDescriptor('wrong-checksum.gbk', 'checksum', {
        checksum: '0'.repeat(64)
      }),
      message: /wrong-checksum \(wrong-checksum\.gbk\) checksum does not match/
    }
  ];

  for (const { id, descriptor, message } of cases) {
    resetMetrics();
    const table = adoptCurrentSessionResources({ [id]: descriptor });
    const view = createSessionResourceFileView(table, id);
    const first = readSessionResourceBytes(view);
    const second = readSessionResourceBytes(view);
    assert.strictEqual(second, first);
    await assert.rejects(first, message);
    await assert.rejects(readFileText(view), message);
    assert.equal(metricTotal('base64DecodeCount', id), 1);
    assert.equal(metricTotal('resourceByteReadCount', id), 1);
  }
});

// PD-OI-121, OV-2619: a file-backed resource whose File has no copy and the browser refuses to read (it changed
// on disk after it was chosen) fails Save with INPUT_UNREADABLE and the approved text, not with
// an "invalid encoded data" error that reads as UNKNOWN.
test('a refused read of a file-backed resource keeps its INPUT_UNREADABLE error', async () => {
  const { diagnosticError, normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  const { fileBackedDescriptor } = await import('../../gbdraw/web/js/services/session-resource-backing.js');
  const refused = diagnosticError('INPUT_UNREADABLE', {}, { stage: 'read' });
  const descriptor = fileBackedDescriptor(
    { name: 'chosen.gb', type: '', size: 4, lastModified: 0 },
    () => Promise.reject(refused)
  );
  const table = adoptCurrentSessionResources({ 'record-1-genbank': descriptor });
  const view = createSessionResourceFileView(table, 'record-1-genbank');
  await assert.rejects(readSessionResourceBytes(view), (error) => error === refused);
  const model = normalizeUserFacingError(refused, { operation: 'session-save' });
  assert.equal(model.code, 'INPUT_UNREADABLE');
  assert.equal(model.summary, 'An input could not be read. Replace or reselect it.');
});

test('ordinary uploaded File objects keep the native read behavior', async () => {
  const file = new File(['native upload'], 'native.txt', { type: 'text/plain' });
  let reads = 0;
  const nativeArrayBuffer = file.arrayBuffer.bind(file);
  file.arrayBuffer = async () => {
    reads += 1;
    return nativeArrayBuffer();
  };

  // Reads under one hold share one read; once the hold ends, the page drops the bytes.
  const releaseReads = holdReadBytes();
  assert.equal(await readFileText(file), 'native upload');
  assert.equal(new TextDecoder().decode(await readFileBytes(file)), 'native upload');
  assert.equal(reads, 1);
  releaseReads();
  // The first read settles once its pieces are hashed; a read after that reads again.
  await readFileContentFacts(file);
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(new TextDecoder().decode(await readFileBytes(file)), 'native upload');
  assert.equal(reads, 2);
});

test('adopted export reuses descriptors and encodes only a replacement', async () => {
  const resources = {
    alpha: encodedDescriptor('alpha.gbk', 'alpha\n'),
    beta: encodedDescriptor('beta.gbk', 'beta\n'),
    unused: encodedDescriptor('unused.txt', 'unused\n')
  };
  const committed = adoptRuntimeCanonicalSession({
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: 'linear',
      records: [
        {
          recordKey: 'alpha-record',
          source: { kind: 'genbank', resourceId: 'alpha' }
        },
        {
          recordKey: 'beta-record',
          source: { kind: 'genbank', resourceId: 'beta' }
        }
      ],
      diagramOptions: {},
      layout: {},
      comparisons: [],
      output: { prefix: 'lazy-export', formats: ['svg'], overwrite: false }
    },
    resources,
    webFiles: {}
  });
  const table = adoptCurrentSessionResources(resources);
  const alpha = createSessionResourceFileView(table, 'alpha');
  const beta = createSessionResourceFileView(table, 'beta');
  const state = {
    files: {},
    linearSeqs: [
      { uid: 'alpha-record', gb: alpha, losat_gencode: 1 },
      { uid: 'beta-record', gb: beta, losat_gencode: 1 }
    ],
    linearComparisonPlan: { edges: [] }
  };

  resetMetrics();
  const unchanged = await buildSessionResources(state, committed, bindingDrawings(state));
  assert.equal(metricTotal('base64DecodeCount'), 0);
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.strictEqual(unchanged.resources.alpha, resources.alpha);
  assert.strictEqual(unchanged.resources.beta, resources.beta);
  // E1 (review m4): a resource that no request or binding names is not
  // written, as the Python writer drops it.
  assert.equal(unchanged.resources.unused, undefined);
  assert.deepEqual(unchanged.resources, { alpha: resources.alpha, beta: resources.beta });

  assert.equal(await readFileText(alpha), 'alpha\n');
  resetMetrics();
  const materializedButUnchanged = await buildSessionResources(state, committed, bindingDrawings(state));
  assert.equal(metricTotal('base64DecodeCount'), 0);
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.strictEqual(materializedButUnchanged.resources.alpha, resources.alpha);

  const replacement = new File(['replacement beta\n'], 'beta.gbk', {
    type: 'text/plain',
    lastModified: 8
  });
  state.linearSeqs[1].gb = replacement;
  resetMetrics();
  const changed = await buildSessionResources(state, committed, bindingDrawings(state));
  const changedId = changed.webFiles.bindings.linearSeqs[1].gb.resourceId;
  assert.equal(metricTotal('base64DecodeCount'), 0);
  // Save streams the replacement's bytes into the Session file; it never makes a base64 text.
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.equal(changed.webFiles.bindings.linearSeqs[0].gb.resourceId, 'alpha');
  assert.notEqual(changedId, 'beta');
  assert.strictEqual(changed.resources.alpha, resources.alpha);
  assert.strictEqual(changed.resources.beta, resources.beta);
  assert.equal(
    Buffer.from(changed.resources[changedId].data, 'base64').toString('utf8'),
    'replacement beta\n'
  );
});

test('export keeps an imported draft file distinct from regenerated bytes under its old ID', async () => {
  const id = 'colors-default-colors-file';
  const original = encodedDescriptor('default-colors.tsv', '# imported colors\nCDS\t#123456\n');
  const regenerated = encodedDescriptor('default-colors.tsv', 'CDS\t#123456\n');
  const resources = { [id]: regenerated };
  const committed = adoptRuntimeCanonicalSession({
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA, mode: 'linear', records: [], comparisons: [],
      diagramOptions: { colors: { defaultColorsFile: { resourceId: id, representation: 'file' } } }
    },
    resources, webFiles: {}
  });
  const file = createSessionResourceFileView(adoptCurrentSessionResources({ [id]: original }), id);
  // The card names its file by ref of the input file store.
  const inputFiles = createInputFileStore();
  const ref = inputFiles.registerFile(file);
  const state = { inputFiles, files: { d_color: ref }, linearSeqs: [], linearComparisonPlan: { edges: [] } };

  const exported = await buildSessionResources(state, committed, bindingDrawings(state));
  const binding = exported.webFiles.bindings.d_color;
  assert.notEqual(binding.resourceId, id);
  assert.strictEqual(exported.renderRequest, committed.renderRequest);
  assert.strictEqual(exported.resources[id], regenerated);
  assert.equal(exported.resources[binding.resourceId].data, original.data);
  assert.notEqual(exported.resources[binding.resourceId].name, regenerated.name);
  assert.equal(binding.name, file.name);
  assert.strictEqual(state.files.d_color, ref);
  assert.deepEqual(resources, { [id]: regenerated });

  const invalid = { ...original, size: original.size + 1 };
  state.files.d_color = inputFiles.registerFile(createSessionResourceFileView(
    adoptCurrentSessionResources({ [id]: invalid }), id
  ));
  const invalidFile = state.files.d_color;
  await assert.rejects(buildSessionResources(state, committed, bindingDrawings(state)), /byte size does not match/);
  assert.strictEqual(state.files.d_color, invalidFile);
  assert.deepEqual(resources, { [id]: regenerated });

  // Conflicting IDs within one adopted table remain invalid.
  const conflicting = adoptRuntimeCanonicalSession({
    ...committed, resources: { [id]: regenerated, [` ${id} `]: original }
  });
  await assert.rejects(buildSessionResources({ files: {} }, conflicting),
    /resource IDs must be unique/);
});

test('Save with one card per file preserves ordered sources across a committed ID collision', async () => {
  const resources = {
    'record-1-genbank': encodedDescriptor('first.gb', 'first\n'),
    second: encodedDescriptor('second.gb', 'second\n')
  };
  const table = adoptCurrentSessionResources(resources);
  const circularSources = Object.keys(resources).map((resourceId, index) => ({
    uid: `record-${index + 1}`, gb: createSessionResourceFileView(table, resourceId), gff: null, fasta: null
  }));
  const committed = adoptRuntimeCanonicalSession({
    renderRequest: { schema: CANONICAL_REQUEST_SCHEMA, records: [{ source: { resourceId: 'record-1-genbank' } }] },
    resources: { 'record-1-genbank': encodedDescriptor('first.gb', 'first\nsecond\n') }
  });
  resetMetrics();
  const saved = await buildSessionResources({ files: {}, circularSources }, committed);
  const bindings = saved.webFiles.bindings.circularSources;
  assert.equal(bindings.length, 2);
  assert.deepEqual(bindings.map(card => card.uid), ['record-1', 'record-2']);
  assert.equal(saved.renderRequest, committed.renderRequest);
  assert.deepEqual(bindings.map(card => saved.resources[card.gb.resourceId].data),
    Object.values(resources).map(part => part.data));
  assert.notEqual(bindings[0].gb.resourceId, 'record-1-genbank');
  assert.equal(Object.keys(saved.resources).length, 3);
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.equal(metricTotal('base64DecodeCount'), 1, 'only the colliding new source needs validation');
});

test('same-length adopted, alternate-encoded and new candidates decode/hash at most once per Save', async () => {
  const resources = { alpha: encodedDescriptor('alpha.gb', 'AA'), beta: encodedDescriptor('beta.gb', 'BB') };
  const alternate = { ...resources.alpha, data: 'Q U E =\n' };
  const table = adoptCurrentSessionResources({ differentId: alternate });
  const file = createSessionResourceFileView(table, 'differentId', { name: 'display.gb' });
  const committed = adoptRuntimeCanonicalSession({ renderRequest: { schema: CANONICAL_REQUEST_SCHEMA, records: [] }, resources });
  resetMetrics();
  const result = await buildSessionResources({
    circularSources: [{ uid: 'record-1', gb: file, gff: null, fasta: null }],
    files: { c_depth: [file, new File(['AA'], 'other.txt'), new File(['CC'], 'new.txt')] }
  }, committed);
  assert.equal(result.webFiles.bindings.circularSources[0].gb.resourceId, 'alpha');
  assert.equal(result.webFiles.bindings.c_depth[0].resourceId, 'alpha');
  assert.equal(result.webFiles.bindings.c_depth[1].resourceId, 'alpha');
  assert.equal(metricTotal('base64DecodeCount'), 3);
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.equal(metricTotal('resourceIdentityHashCount'), 5);
  for (const id of ['alpha', 'beta', 'differentId']) assert.equal(metricTotal('resourceByteReadCount', id), 1);
});

// R1-M3 (R6): the checks of a resource descriptor and of a Web input inventory
// refuse with the typed Session-fields error, not a plain Error with their text.
test('a resource descriptor or Web binding the check refuses is the typed Session-fields error (R1-M3)', () => {
  const sessionFields = (error) => error?.code === 'INPUT_INVALID' && error.context?.field === 'schema'
    && error.context?.reason === 'FIELDS';
  assert.throws(() => adoptCurrentSessionResources({ alpha: 5 }), sessionFields);
  assert.throws(() => validateWebFileBindings({ bindings: { schema: 2, c_gb: { resourceId: 'missing' } } }, {}), sessionFields);
});
