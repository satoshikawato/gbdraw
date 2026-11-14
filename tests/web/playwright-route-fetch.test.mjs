// OV-713 (2026-10-13): the local static server can reset a connection under CI
// load, and a browser test that served the request through route.fetch() then
// failed with "route.fetch: read ECONNRESET". Browser tests fetch routed
// requests through fetchRoute (tests/web/helpers/app-lifecycle.cjs), which
// retries once on a connection error and rethrows any other error unchanged.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import test from 'node:test';

const { fetchRoute, routeCrossOriginIsolated } = createRequire(import.meta.url)('./helpers/app-lifecycle.cjs');

const fakeRoute = (...outcomes) => {
  const calls = [];
  return {
    calls,
    fetch: async (options) => {
      calls.push(options);
      const outcome = outcomes[calls.length - 1];
      if (outcome instanceof Error) throw outcome;
      return outcome;
    }
  };
};

test('fetchRoute retries once on a connection error', async () => {
  for (const message of [
    'route.fetch: read ECONNRESET',
    'route.fetch: connect ECONNREFUSED 127.0.0.1:4173',
    'route.fetch: socket hang up'
  ]) {
    const route = fakeRoute(new Error(message), 'response');
    assert.equal(await fetchRoute(route, { timeout: 5 }), 'response', message);
    assert.deepEqual(route.calls, [{ timeout: 5 }, { timeout: 5 }], message);
  }
});

test('fetchRoute rethrows a second connection error and any other error unchanged', async () => {
  const second = new Error('route.fetch: read ECONNRESET');
  const twice = fakeRoute(new Error('route.fetch: read ECONNRESET'), second, 'unused');
  await assert.rejects(fetchRoute(twice), (error) => error === second);
  assert.equal(twice.calls.length, 2);

  const other = new Error('route.fetch: Test ended');
  const once = fakeRoute(other, 'unused');
  await assert.rejects(fetchRoute(once), (error) => error === other);
  assert.equal(once.calls.length, 1);
});

// Run 43572887428 (2026-11-14): a font request still in flight when its test
// ended failed the shard outside any test with "route.fetch: Test ended".
test('the cross-origin isolated route drops a request that outlives its test and rethrows other errors', async () => {
  let serve;
  await routeCrossOriginIsolated({ route: async (_pattern, handler) => { serve = handler; } });
  const routeFailingWith = (message) => ({
    ...fakeRoute(new Error(message)),
    request: () => ({ url: () => 'http://127.0.0.1:4173/gbdraw/web/vendor/fonts/font.woff2' }),
    fulfill: async () => assert.fail('fulfilled a failed fetch')
  });
  for (const message of [
    'route.fetch: Test ended.',
    'route.fetch: Target page, context or browser has been closed'
  ]) {
    assert.equal(await serve(routeFailingWith(message)), undefined, message);
  }
  await assert.rejects(serve(routeFailingWith('route.fetch: net::ERR_FAILED')), /ERR_FAILED/);
});

test('browser tests call route.fetch only through fetchRoute', () => {
  const owner = 'tests/web/helpers/app-lifecycle.cjs';
  const files = execFileSync('git', ['ls-files', '-z', '--', 'tests/web'], { encoding: 'utf8' })
    .split('\0')
    .filter((path) => /\.(?:c|m)?js$/.test(path) && !path.endsWith('.test.mjs') && path !== owner);
  const offenders = files.filter((path) => /\broute\.fetch\(/.test(readFileSync(path, 'utf8')));
  assert.deepEqual(offenders, [], `Use fetchRoute from ${owner} instead of route.fetch():\n${offenders.join('\n')}`);
});
