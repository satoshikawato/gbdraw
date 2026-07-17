import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';

const workerSource = await readFile(
  join(process.cwd(), 'gbdraw', 'web', 'cloudflare-worker.js'),
  'utf8'
);
const workerModuleUrl = `data:text/javascript;base64,${Buffer.from(workerSource).toString('base64')}`;
const { default: worker } = await import(workerModuleUrl);

const assetRequests = [];
const env = {
  ASSETS: {
    fetch: async (request) => {
      const url = new URL(request.url);
      assetRequests.push({ method: request.method, pathname: url.pathname, search: url.search });
      return new Response(`asset:${url.pathname}`, {
        status: 200,
        headers: { 'content-type': 'text/plain' },
      });
    },
  },
};

const dispatch = (path, options = {}) =>
  worker.fetch(new Request(`https://gbdraw.test${path}`, options), env);

{
  const response = await dispatch('/gallery/palettes?palette=orchid');
  assert.equal(response.status, 308);
  assert.equal(response.headers.get('location'), 'https://gbdraw.test/gallery/palettes/?palette=orchid');
  assert.equal(assetRequests.length, 0);
}

{
  const response = await dispatch('/gallery/palettes/?palette=orchid');
  assert.equal(response.status, 200);
  assert.equal(await response.text(), 'asset:/gallery/palettes/');
  assert.deepEqual(assetRequests.pop(), {
    method: 'GET',
    pathname: '/gallery/palettes/',
    search: '',
  });
}

{
  const response = await dispatch('/gallery/HmmtDNA_basic_circular#Tutorial');
  assert.equal(response.status, 200);
  assert.equal(await response.text(), 'asset:/gallery/');
  assert.deepEqual(assetRequests.pop(), {
    method: 'GET',
    pathname: '/gallery/',
    search: '',
  });
}

{
  const response = await dispatch('/gallery/HmmtDNA_basic_circular/?preview=1');
  assert.equal(response.status, 308);
  assert.equal(
    response.headers.get('location'),
    'https://gbdraw.test/gallery/HmmtDNA_basic_circular?preview=1'
  );
  assert.equal(assetRequests.length, 0);
}

for (const path of [
  '/gallery/examples',
  '/gallery/palettes/palette-explorer.js',
  '/gallery/tutorials/HmmtDNA_basic_circular.json',
]) {
  const response = await dispatch(path);
  assert.equal(response.status, 200);
  assert.equal(await response.text(), `asset:${path}`);
  assert.equal(assetRequests.pop().pathname, path);
}

{
  const response = await dispatch('/gallery/HmmtDNA_basic_circular', { method: 'POST' });
  assert.equal(response.status, 200);
  assert.equal(await response.text(), 'asset:/gallery/HmmtDNA_basic_circular');
  assert.equal(assetRequests.pop().pathname, '/gallery/HmmtDNA_basic_circular');
}

{
  const response = await dispatch('/gallery/palettes/', { method: 'HEAD' });
  assert.equal(response.status, 200);
  assert.deepEqual(assetRequests.pop(), {
    method: 'HEAD',
    pathname: '/gallery/palettes/',
    search: '',
  });
}

assert.equal(assetRequests.length, 0);
