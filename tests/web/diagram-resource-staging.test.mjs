import assert from 'node:assert/strict';
import test from 'node:test';

const metrics = [];
const lifecycle = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric: (metric) => metrics.push(metric),
  onSessionLifecycleEvent: (event) => lifecycle.push(event)
};

const { createDiagramResourceTransport } = await import(
  '../../gbdraw/web/js/services/diagram-resource-staging.js'
);

const descriptor = (text, name = 'record.gb') => ({
  kind: 'genbank',
  name,
  type: 'text/plain',
  size: Buffer.byteLength(text),
  lastModified: 1,
  encoding: 'base64',
  data: Buffer.from(text).toString('base64')
});

test('render staging selects references, transfers bytes once, and reuses the cache', async () => {
  const transport = createDiagramResourceTransport();
  const record = descriptor('LOCUS record\nORIGIN\n//\n');
  const unused = descriptor('unused', 'unused.txt');
  const request = {
    records: [{ source: { kind: 'genbank', resourceId: 'record-1-genbank' } }],
    diagramOptions: {}
  };
  const resources = {
    'record-1-genbank': record,
    'unused-resource': unused
  };

  const first = await transport.prepare({ request, resources });
  assert.deepEqual(
    first.resourceManifest.map(({ resourceId }) => resourceId),
    ['record-1-genbank']
  );
  assert.equal(first.stagedResources.length, 1);
  assert.equal(
    Buffer.from(first.stagedResources[0].bytes).toString('utf8'),
    'LOCUS record\nORIGIN\n//\n'
  );
  first.commit();

  const second = await transport.prepare({ request, resources });
  assert.equal(second.stagedResources.length, 0);
  assert.equal(
    second.resourceManifest[0].cacheToken,
    first.resourceManifest[0].cacheToken
  );
  second.commit();

  const replacement = descriptor('LOCUS other!\nORIGIN\n//\n');
  const third = await transport.prepare({
    request,
    resources: { ...resources, 'record-1-genbank': replacement }
  });
  assert.equal(third.stagedResources.length, 1);
  assert.notEqual(
    third.resourceManifest[0].cacheToken,
    first.resourceManifest[0].cacheToken
  );

  assert.equal(
    metrics.filter(({ name }) => name === 'resourceMaterializationCount').length,
    2
  );
  assert.equal(
    metrics.filter(({ name }) => name === 'workerResourceCacheHitCount').length,
    1
  );
  assert.equal(
    lifecycle.filter(({ name }) => name === 'resource-stage-end').length,
    3
  );
});
