import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { createHash, webcrypto } from 'node:crypto';
import test from 'node:test';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric(metric) {
    metrics.push({ ...metric });
  }
};

const {
  adoptCurrentSessionResources,
  createCombinedSessionResourceFileView,
  createSessionResourceFileView,
  readSessionResourceBytes,
  readSessionResourceText
} = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const {
  readFileBytes,
  readFileText
} = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { adoptRuntimeCanonicalSession } = await import(
  '../../gbdraw/web/js/services/session-authority.js'
);
const { buildSessionResources } = await import(
  '../../gbdraw/web/js/services/session-resources.js'
);

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

test('combined file views report exact canonical byte sizes without preview decoding', async () => {
  const cases = [
    {
      name: 'both resources end with LF',
      parts: ['a\n', 'ab\n'],
      expected: 'a\nab\n'
    },
    {
      name: 'neither resource ends with LF',
      parts: ['a', 'abc'],
      expected: 'a\nabc\n'
    },
    {
      name: 'one resource ends with LF',
      parts: ['alpha\n', 'be'],
      expected: 'alpha\nbe\n'
    }
  ];

  for (const { name, parts, expected } of cases) {
    resetMetrics();
    const resources = Object.fromEntries(parts.map((text, index) => [
      `part-${index + 1}`,
      encodedDescriptor(`part-${index + 1}.gbk`, text)
    ]));
    const table = adoptCurrentSessionResources(resources);
    const view = createCombinedSessionResourceFileView(table, Object.keys(resources));
    const expectedBytes = new TextEncoder().encode(expected);

    assert.equal(metrics.length, 0, `${name}: metadata projection stays lazy`);
    assert.equal(view.size, expectedBytes.byteLength, `${name}: declared size is exact`);

    const backingRead = readSessionResourceBytes(view);
    assert.strictEqual(
      readSessionResourceBytes(view),
      backingRead,
      `${name}: concurrent combined backing reads share one Promise`
    );
    const genericRead = readFileBytes(view);
    assert.strictEqual(
      genericRead,
      backingRead,
      `${name}: generic and combined boundaries share one Promise`
    );
    assert.strictEqual(readFileBytes(view), genericRead);

    const [backingBytes, genericBytes] = await Promise.all([backingRead, genericRead]);
    assert.strictEqual(
      genericBytes,
      backingBytes,
      `${name}: generic and combined boundaries share one Uint8Array`
    );
    assert.equal(view.size, genericBytes.byteLength, `${name}: materialized size is exact`);
    assert.equal(new TextDecoder().decode(genericBytes), expected);

    Object.keys(resources).forEach((resourceId) => {
      assert.equal(metricTotal('base64DecodeCount', resourceId), 1);
      assert.equal(metricTotal('resourceByteReadCount', resourceId), 1);
      assert.equal(
        metricTotal('decodedByteCount', resourceId),
        resources[resourceId].size
      );
    });
  }
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

test('ordinary uploaded File objects keep the native read behavior', async () => {
  const file = new File(['native upload'], 'native.txt', { type: 'text/plain' });
  let reads = 0;
  const nativeArrayBuffer = file.arrayBuffer.bind(file);
  file.arrayBuffer = async () => {
    reads += 1;
    return nativeArrayBuffer();
  };

  assert.equal(await readFileText(file), 'native upload');
  assert.equal(new TextDecoder().decode(await readFileBytes(file)), 'native upload');
  assert.equal(reads, 1);
});

test('adopted export reuses descriptors and encodes only a replacement', async () => {
  const resources = {
    alpha: encodedDescriptor('alpha.gbk', 'alpha\n'),
    beta: encodedDescriptor('beta.gbk', 'beta\n'),
    unused: encodedDescriptor('unused.txt', 'unused\n')
  };
  const committed = adoptRuntimeCanonicalSession({
    renderRequest: {
      schema: 5,
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
  const unchanged = await buildSessionResources(state, committed);
  assert.equal(metricTotal('base64DecodeCount'), 0);
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.strictEqual(unchanged.resources.alpha, resources.alpha);
  assert.strictEqual(unchanged.resources.beta, resources.beta);
  assert.strictEqual(unchanged.resources.unused, resources.unused);
  assert.deepEqual(unchanged.resources, resources);

  assert.equal(await readFileText(alpha), 'alpha\n');
  resetMetrics();
  const materializedButUnchanged = await buildSessionResources(state, committed);
  assert.equal(metricTotal('base64DecodeCount'), 0);
  assert.equal(metricTotal('base64EncodeCount'), 0);
  assert.strictEqual(materializedButUnchanged.resources.alpha, resources.alpha);

  const replacement = new File(['replacement beta\n'], 'beta.gbk', {
    type: 'text/plain',
    lastModified: 8
  });
  state.linearSeqs[1].gb = replacement;
  resetMetrics();
  const changed = await buildSessionResources(state, committed);
  const changedId = changed.webFiles.bindings.linearSeqs[1].gb.resourceId;
  assert.equal(metricTotal('base64DecodeCount'), 0);
  assert.equal(metricTotal('base64EncodeCount'), 1);
  assert.equal(metricTotal('encodedByteCount'), replacement.size);
  assert.equal(changed.webFiles.bindings.linearSeqs[0].gb.resourceId, 'alpha');
  assert.notEqual(changedId, 'beta');
  assert.strictEqual(changed.resources.alpha, resources.alpha);
  assert.strictEqual(changed.resources.beta, resources.beta);
  assert.equal(
    Buffer.from(changed.resources[changedId].data, 'base64').toString('utf8'),
    'replacement beta\n'
  );
});
