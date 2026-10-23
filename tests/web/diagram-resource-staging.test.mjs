import assert from 'node:assert/strict';
import test from 'node:test';
import { createHash } from 'node:crypto';

const metrics = [];
const lifecycle = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric: (metric) => metrics.push(metric),
  onSessionLifecycleEvent: (event) => lifecycle.push(event)
};

const { createDiagramResourceTransport } = await import(
  '../../gbdraw/web/js/services/diagram-resource-staging.js'
);
const { setResourcePayloadOwner } = await import(
  '../../gbdraw/web/js/services/resource-payload-owner.js'
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

  const fileBackedRecord = setResourcePayloadOwner(
    { ...record },
    { arrayBuffer: async () => Buffer.from('LOCUS record\nORIGIN\n//\n').buffer }
  );
  const second = await transport.prepare({
    request,
    resources: { ...resources, 'record-1-genbank': fileBackedRecord }
  });
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

test('unacknowledged preparation, replacement, and reset never publish stale ownership', async () => {
  const transport = createDiagramResourceTransport();
  const request = { records: [{ source: { kind: 'genbank', resourceId: 'source' } }] };
  const prepare = text => transport.prepare({ request, resources: { source: descriptor(text) } });
  const abandoned = await prepare('A');
  const first = await prepare('A');
  assert.equal(first.stagedResources.length, 1, 'prepare alone does not acknowledge bytes');
  assert.equal(first.commit(), true);
  assert.equal(abandoned.commit(), false, 'older preparation cannot replace an acknowledged state');
  assert.equal(first.commit(), false);
  const pending = await prepare('B');
  assert.equal((await prepare('A')).stagedResources.length, 0);
  assert.equal((await prepare('B')).stagedResources.length, 1);
  assert.equal(pending.commit(), true);
  assert.equal((await prepare('A')).stagedResources.length, 1);
  assert.equal((await prepare('B')).stagedResources.length, 0);
  const beforeReset = await prepare('B');
  transport.reset();
  assert.equal(beforeReset.commit(), false, 'reset invalidates outstanding preparations');
  assert.equal((await prepare('B')).stagedResources.length, 1);
});

const genbankRequest = (...resourceIds) => ({
  records: resourceIds.map((resourceId) => ({ source: { kind: 'genbank', resourceId } }))
});
const sha256Hex = (text) => createHash('sha256').update(text).digest('hex');
const stagedTokens = (prepared) => prepared.stagedResources.map(({ cacheToken }) => cacheToken);

// A4 (IN-A): the Worker keeps the resources of every drawing's committed request, named by the sha256 of their bytes.
test('alternating drawings reuse the resources every drawing committed', async () => {
  const committed = [];
  const transport = createDiagramResourceTransport({ readCommittedResourceSets: () => committed });
  const circular = { request: genbankRequest('circular-genbank'), resources: { 'circular-genbank': descriptor('CIRCULAR') } };
  const linear = {
    request: genbankRequest('linear-a', 'linear-b'),
    resources: { 'linear-a': descriptor('LINEAR A'), 'linear-b': descriptor('LINEAR B') }
  };

  const first = await transport.prepare(circular);
  assert.deepEqual(stagedTokens(first), [sha256Hex('CIRCULAR')]);
  assert.equal(first.commit(), true);
  committed.push(circular);

  const second = await transport.prepare(linear);
  assert.deepEqual(stagedTokens(second), [sha256Hex('LINEAR A'), sha256Hex('LINEAR B')]);
  assert.deepEqual(second.retainedResourceTokens, [sha256Hex('CIRCULAR')], 'the Circular drawing keeps its resource');
  assert.equal(second.commit(), true);
  committed.push(linear);

  const third = await transport.prepare(circular);
  assert.deepEqual(third.stagedResources, [], 'Circular again transfers nothing');
  assert.deepEqual(third.resourceManifest.map(({ cacheToken }) => cacheToken), [sha256Hex('CIRCULAR')]);
  assert.deepEqual([...third.retainedResourceTokens].sort(), [sha256Hex('LINEAR A'), sha256Hex('LINEAR B')].sort());
  assert.equal(third.commit(), true);

  const fourth = await transport.prepare(linear);
  assert.deepEqual(fourth.stagedResources, [], 'Linear again transfers nothing');
});

test('equal bytes under two resource IDs share one token and one transfer', async () => {
  const transport = createDiagramResourceTransport();
  const text = 'LOCUS same\nORIGIN\n//\n';
  const fileBacked = setResourcePayloadOwner(
    descriptor(text, 'copy.gb'),
    { arrayBuffer: async () => Buffer.from(text).buffer }
  );
  const prepared = await transport.prepare({
    request: genbankRequest('first-copy', 'second-copy'),
    resources: { 'first-copy': descriptor(text), 'second-copy': fileBacked }
  });
  assert.deepEqual(
    prepared.resourceManifest.map(({ resourceId, cacheToken }) => [resourceId, cacheToken]),
    [['first-copy', sha256Hex(text)], ['second-copy', sha256Hex(text)]]
  );
  assert.deepEqual(stagedTokens(prepared), [sha256Hex(text)]);
});

test('a resource no committed drawing names is released', async () => {
  const committed = [];
  const transport = createDiagramResourceTransport({ readCommittedResourceSets: () => committed });
  const circular = { request: genbankRequest('circular-genbank'), resources: { 'circular-genbank': descriptor('OLD CIRCULAR') } };
  const linear = { request: genbankRequest('linear-genbank'), resources: { 'linear-genbank': descriptor('LINEAR') } };
  (await transport.prepare(circular)).commit();
  committed.push(circular);
  (await transport.prepare(linear)).commit();
  committed.push(linear);
  // The Circular drawing now commits other bytes: nothing names the old ones.
  committed[0] = { request: circular.request, resources: { 'circular-genbank': descriptor('NEW CIRCULAR') } };
  const next = await transport.prepare(linear);
  assert.deepEqual(next.retainedResourceTokens, [], 'the old Circular bytes are not retained');
  assert.equal(next.commit(), true);
  assert.deepEqual(stagedTokens(await transport.prepare(circular)), [sha256Hex('OLD CIRCULAR')]);
});

test('staging decodes base64 with the one shared decoder (OV-1560)', async () => {
  const bytes = Uint8Array.from({ length: 200_003 }, (_, index) => (index * 89 + 7) % 256);
  const resources = {
    big: { kind: 'genbank', name: 'big.gb', type: '', size: bytes.length, lastModified: 0,
      encoding: 'base64', data: Buffer.from(bytes).toString('base64') },
    bad: { kind: 'genbank', name: 'bad.gb', type: '', size: 3, lastModified: 0, encoding: 'base64', data: 'TQ!A' }
  };
  const request = (resourceId) => ({ records: [{ source: { kind: 'genbank', resourceId } }] });
  const native = Object.getOwnPropertyDescriptor(Uint8Array.prototype, 'setFromBase64');
  const decodeWith = async (setFromBase64) => {
    Object.defineProperty(Uint8Array.prototype, 'setFromBase64', { ...native, value: setFromBase64 });
    try {
      return await createDiagramResourceTransport().prepare({ request: request('big'), resources });
    } finally {
      Object.defineProperty(Uint8Array.prototype, 'setFromBase64', native);
    }
  };
  let nativeCalls = 0;
  const viaNative = await decodeWith(function (...args) { nativeCalls += 1; return native.value.apply(this, args); });
  assert.ok(nativeCalls > 0, 'the engine decoder writes the transfer buffer');
  assert.deepEqual(new Uint8Array(viaNative.stagedResources[0].bytes), bytes);
  const viaFallback = await decodeWith(undefined);
  assert.deepEqual(new Uint8Array(viaFallback.stagedResources[0].bytes), bytes);
  await assert.rejects(
    createDiagramResourceTransport().prepare({ request: request('bad'), resources }),
    /Canonical resource bad contains invalid base64 data\./
  );
});

test('two Linear cards that read one Session resource both find its committed resource id (OV-1555)', async () => {
  const { resourceIdOwnedBy } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
  const { adoptCurrentSessionResources, createSessionResourceFileView, sessionResourceSource } = await import(
    '../../gbdraw/web/js/services/session-resource-backing.js'
  );
  const table = adoptCurrentSessionResources({ 'record-1-genbank': descriptor('LOCUS two\nORIGIN\n//\n') });
  const firstCard = createSessionResourceFileView(table, 'record-1-genbank');
  const secondCard = createSessionResourceFileView(table, 'record-1-genbank');
  // A Save serializes each card's view to the one source descriptor, and the last card names its owner;
  // the request keeps one resource for both cards (`session-request.js` `addFile`).
  const source = sessionResourceSource(firstCard).descriptor;
  setResourcePayloadOwner(source, firstCard);
  setResourcePayloadOwner(source, secondCard);
  const committed = {
    'record-1-genbank': setResourcePayloadOwner({ ...source }, secondCard)
  };
  assert.equal(resourceIdOwnedBy(committed, firstCard, 'genbank'), 'record-1-genbank');
  assert.equal(resourceIdOwnedBy(committed, secondCard, 'genbank'), 'record-1-genbank');
  assert.equal(resourceIdOwnedBy(committed, firstCard, 'fasta'), null);
  const uploaded = { name: 'new.gb', size: 3 };
  assert.equal(resourceIdOwnedBy(committed, uploaded, 'genbank'), null);
  const uploadedResources = { 'record-2-genbank': setResourcePayloadOwner(descriptor('LOCUS new'), uploaded) };
  assert.equal(resourceIdOwnedBy(uploadedResources, uploaded, 'genbank'), 'record-2-genbank');
});
