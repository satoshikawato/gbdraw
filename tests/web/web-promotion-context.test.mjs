import assert from 'node:assert/strict';
import test from 'node:test';

import { classifyWebChangeContext } from '../../tools/web-promotion-context.mjs';

const BASE_SHA = '1'.repeat(40);
const HEAD_SHA = '2'.repeat(40);
const REPOSITORY = 'satoshikawato/gbdraw';

const pullRequestPayload = ({
  baseRef = 'main',
  headRef = 'dev',
  baseSha = BASE_SHA,
  headSha = HEAD_SHA,
  baseRepository = REPOSITORY,
  headRepository = REPOSITORY
} = {}) => JSON.stringify({
  pull_request: {
    base: {
      ref: baseRef,
      sha: baseSha,
      repo: { full_name: baseRepository }
    },
    head: {
      ref: headRef,
      sha: headSha,
      repo: { full_name: headRepository }
    }
  }
});

const classify = (overrides = {}) => classifyWebChangeContext({
  eventName: 'pull_request',
  currentRepository: REPOSITORY,
  eventPayloadSource: pullRequestPayload(),
  baseSha: BASE_SHA,
  headSha: HEAD_SHA,
  ...overrides
});

test('same-repository dev to main pull requests are promotions', () => {
  for (const eventName of ['pull_request', 'pull_request_target']) {
    const result = classify({ eventName });
    assert.equal(result.context, 'PROMOTION');
    assert.equal(result.isPromotion, true);
    assert.deepEqual(result.errors, []);
  }
});

test('fork dev to main pull requests are invalid promotions', () => {
  const result = classify({
    eventPayloadSource: pullRequestPayload({ headRepository: 'contributor/gbdraw' })
  });
  assert.equal(result.context, 'PROMOTION');
  assert.equal(result.isPromotion, true);
  assert.deepEqual(
    result.errors,
    ['Promotion pull requests must use dev from the current repository.']
  );
});

test('same-repository feature to main pull requests are invalid promotions', () => {
  const result = classify({
    eventPayloadSource: pullRequestPayload({ headRef: 'feature/example' })
  });
  assert.equal(result.context, 'PROMOTION');
  assert.equal(result.isPromotion, true);
  assert.deepEqual(
    result.errors,
    ['Promotion pull requests must use dev as the head branch.']
  );
});

test('same-repository dev to release pull requests are ordinary changes', () => {
  const result = classify({
    eventPayloadSource: pullRequestPayload({ baseRef: 'release' })
  });
  assert.equal(result.context, 'ORDINARY');
  assert.equal(result.isPromotion, false);
  assert.deepEqual(result.errors, []);
});

test('checker SHAs must match the pull request payload', () => {
  const baseMismatch = classify({ baseSha: '3'.repeat(40) });
  assert.equal(baseMismatch.context, 'PROMOTION');
  assert.equal(baseMismatch.isPromotion, true);
  assert.deepEqual(
    baseMismatch.errors,
    ['Checker base SHA does not match the GitHub event payload.']
  );

  const headMismatch = classify({ headSha: '4'.repeat(40) });
  assert.equal(headMismatch.context, 'PROMOTION');
  assert.equal(headMismatch.isPromotion, true);
  assert.deepEqual(
    headMismatch.errors,
    ['Checker head SHA does not match the GitHub event payload.']
  );
});

test('missing or malformed pull request payloads fail closed', () => {
  for (const eventPayloadSource of ['', '{', '{}', '[]']) {
    const result = classify({ eventPayloadSource });
    assert.equal(result.context, 'ORDINARY');
    assert.equal(result.isPromotion, false);
    assert.ok(result.errors.length, eventPayloadSource);
  }
});

test('an arbitrary promotion variable cannot select promotion context', () => {
  const result = classifyWebChangeContext({
    eventName: 'push',
    currentRepository: REPOSITORY,
    eventPayloadSource: pullRequestPayload(),
    baseSha: BASE_SHA,
    headSha: HEAD_SHA,
    WEB_PROMOTION: 'true'
  });
  assert.equal(result.context, 'ORDINARY');
  assert.equal(result.isPromotion, false);
  assert.deepEqual(result.errors, []);
});
