const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream'
};

let server;
let baseUrl;

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
      if ((!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) || !existsSync(filePath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(filePath)] || 'application/octet-stream'
      });
      createReadStream(filePath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      baseUrl = `http://127.0.0.1:${server.address().port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('Linear record rows and N-to-M comparison batches remain keyed by sequence uid', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const state = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    app.addLinearSeq();
    app.addLinearSeq();
    app.linearRecordLayoutEnabled = true;
    app.syncLinearRecordLayout();
    app.setLinearRecordRow(app.linearSeqs[0].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[1].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[2].uid, 2);
    app.setLinearRecordRow(app.linearSeqs[3].uid, 2);
    app.addLinearComparisonBatch(true);
    return {
      tokens: app.linearLayoutTokens,
      comparisonCount: app.linearComparisons.length,
      endpoints: app.linearComparisons.map((item) => [item.queryUid, item.subjectUid]),
      uniqueUids: new Set(app.linearSeqs.map((item) => item.uid)).size
    };
  });

  expect(state.tokens).toEqual(['#1@1', '#2@1', '#3@2', '#4@2']);
  expect(state.comparisonCount).toBe(4);
  expect(state.uniqueUids).toBe(4);
  expect(new Set(state.endpoints.map((pair) => pair.join('->'))).size).toBe(4);
  await expect(page.locator('input[aria-label="Linear record row"]')).toHaveCount(4);
  await expect(page.getByText('All adjacent-row pairs', { exact: true })).toBeVisible();
});
