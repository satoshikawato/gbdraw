// T5 (2026-10-02): Chromium's inspector can collect a page promise that
// page.evaluate awaits for a long time; Playwright then reports "Execution
// context was destroyed" although no navigation happened. Browser tests start
// long app operations through evaluateWithRetainedPromise
// (tests/web/helpers/app-lifecycle.cjs), which starts the operation, polls for
// its settlement, and rethrows a rejection. This guard keeps direct awaits out.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import test from 'node:test';
import vm from 'node:vm';

const { evaluateWithRetainedPromise } = createRequire(import.meta.url)('./helpers/app-lifecycle.cjs');

// App operations that run a diagram render or a whole Session import or save.
const LONG_OPERATIONS = ['runAnalysis', 'importSession', 'saveSessionWithTitle'];
const OPERATION = `[\\w$.?]*(?:${LONG_OPERATIONS.join('|')})\\(`;
const AWAITED_OPERATION = [
  new RegExp(`\\b(?:await|return)\\s+${OPERATION}`),
  new RegExp(`=>\\s*\\(?\\s*${OPERATION}`)
];
const STORED_OPERATION = new RegExp(`window\\.([\\w$]+)\\s*=\\s*${OPERATION}`, 'g');

const skipQuoted = (source, index) => {
  const quote = source[index];
  let cursor = index + 1;
  while (cursor < source.length && source[cursor] !== quote) {
    if (source[cursor] === '\\') cursor += 1;
    else if (quote === '`' && source.startsWith('${', cursor)) cursor = skipBalanced(source, cursor + 1) - 1;
    cursor += 1;
  }
  return cursor + 1;
};

const skipRegex = (source, index) => {
  let cursor = index + 1;
  let inClass = false;
  while (cursor < source.length && source[cursor] !== '\n') {
    const char = source[cursor];
    if (char === '\\') cursor += 1;
    else if (char === '[') inClass = true;
    else if (char === ']') inClass = false;
    else if (char === '/' && !inClass) break;
    cursor += 1;
  }
  cursor += 1;
  while (/[a-z]/i.test(source[cursor] || '')) cursor += 1;
  return cursor;
};
const REGEX_PREFIX = /(?:^|[(,=:[!&|?{};+\-*%<>~^]|\breturn|\btypeof|\bcase)\s*$/;

// Returns the index after the bracket that closes the one at `start`, skipping
// strings, template literals, comments, and regular expression literals.
function skipBalanced(source, start) {
  let depth = 0;
  let cursor = start;
  while (cursor < source.length) {
    const char = source[cursor];
    if (char === "'" || char === '"' || char === '`') {
      cursor = skipQuoted(source, cursor);
    } else if (source.startsWith('//', cursor)) {
      const lineEnd = source.indexOf('\n', cursor);
      cursor = lineEnd < 0 ? source.length : lineEnd;
    } else if (source.startsWith('/*', cursor)) {
      const commentEnd = source.indexOf('*/', cursor + 2);
      cursor = commentEnd < 0 ? source.length : commentEnd + 2;
    } else if (char === '/' && REGEX_PREFIX.test(source.slice(Math.max(0, cursor - 16), cursor))) {
      cursor = skipRegex(source, cursor);
    } else {
      if ('([{'.includes(char)) depth += 1;
      if (')]}'.includes(char)) {
        depth -= 1;
        if (depth === 0) return cursor + 1;
      }
      cursor += 1;
    }
  }
  return source.length;
}

const findAwaitedLongOperations = (source) => {
  const stored = [...source.matchAll(STORED_OPERATION)].map((match) => match[1]);
  const storedAwait = stored.length
    ? new RegExp(`\\b(?:await|return)\\s+window\\.(?:${stored.map((name) => name.replace(/\$/g, '\\$')).join('|')})\\b`)
    : null;
  const findings = [];
  for (const match of source.matchAll(/\.(?:evaluate|evaluateHandle)\(/g)) {
    const open = match.index + match[0].length - 1;
    const call = source.slice(open, skipBalanced(source, open));
    if (AWAITED_OPERATION.some((pattern) => pattern.test(call)) || storedAwait?.test(call)) {
      findings.push(source.slice(0, match.index).split('\n').length);
    }
  }
  return findings;
};

test('the detector flags awaited long app operations and accepts retained or started ones', () => {
  const flagged = [
    'await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());',
    'await page.evaluate(async () => { const result = await app.runAnalysis(); return result; });',
    'await page.evaluate(async () => ({ result: await window.__GBDRAW_APP__.runAnalysis() }));',
    'await page.evaluate((file) => { return window.__GBDRAW_APP__.importSession(file); }, file);',
    'await run.page.evaluate(async () => { await window.__GBDRAW_APP__.saveSessionWithTitle(); });',
    [
      'await page.evaluate(() => { window.__RUN__ = window.__GBDRAW_APP__.runAnalysis(); });',
      'await page.evaluate(async () => ({ result: await window.__RUN__ }));'
    ].join('\n')
  ];
  for (const source of flagged) {
    assert.equal(findAwaitedLongOperations(source).length, 1, source);
  }
  const accepted = [
    'await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());',
    'await evaluateWithRetainedPromise(page, async () => ({ result: await window.__RUN__ }));',
    'await page.evaluate(() => { window.__RUN__ = window.__GBDRAW_APP__.runAnalysis(); });',
    'await page.evaluate(() => { window.__GBDRAW_APP__.runAnalysis().then(done, fail); });',
    'await page.evaluate(() => window.__GBDRAW_APP__.results.length); // await runAnalysis()',
    [
      "await page.evaluate(async ({ path }) => fetch(`/${path.replace(/^\\.\\//, '')}`), example);",
      'await evaluateWithRetainedPromise(page, async () => { await app.runAnalysis(); });'
    ].join('\n')
  ];
  for (const source of accepted) {
    assert.deepEqual(findAwaitedLongOperations(source), [], source);
  }
});

test('browser tests do not await long app operations inside page.evaluate', () => {
  const files = execFileSync('git', ['ls-files', '-z', '--', 'tests/web'], { encoding: 'utf8' })
    .split('\0')
    .filter((path) => /\.(?:c|m)?js$/.test(path) && !path.endsWith('.test.mjs'));
  assert.ok(files.length > 100, `expected the browser test inventory, found ${files.length} files`);
  const offenders = files.flatMap((path) => (
    findAwaitedLongOperations(readFileSync(path, 'utf8')).map((line) => `${path}:${line}`)
  ));
  assert.deepEqual(
    offenders,
    [],
    `Start ${LONG_OPERATIONS.join(', ')} with evaluateWithRetainedPromise from `
      + 'tests/web/helpers/app-lifecycle.cjs instead of awaiting it inside page.evaluate:\n'
      + offenders.join('\n')
  );
});

// A page stand-in that runs each evaluation in one shared window context.
const createPage = (app) => {
  const context = vm.createContext({ Promise, Error });
  context.window = context;
  context.__GBDRAW_APP__ = app;
  const run = (callback, argument) => vm.runInContext(
    typeof callback === 'string'
      ? callback
      : `(${callback.toString()})(${JSON.stringify(argument) ?? 'undefined'})`,
    context
  );
  return {
    context,
    isClosed: () => false,
    evaluate: async (callback, argument) => run(callback, argument),
    waitForFunction: async (callback, argument) => {
      while (!run(callback, argument)) await new Promise((resolve) => setTimeout(resolve, 1));
    }
  };
};

test('retained evaluation starts the operation once, returns its value, and rethrows its rejection', async () => {
  let starts = 0;
  let settle;
  const page = createPage({
    runAnalysis: () => {
      starts += 1;
      return new Promise((resolve, reject) => { settle = { resolve, reject }; });
    }
  });
  const resolved = evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  await new Promise((resolve) => setTimeout(resolve, 5));
  settle.resolve({ status: 'ok' });
  assert.deepEqual(await resolved, { status: 'ok' });

  const rejected = evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  await new Promise((resolve) => setTimeout(resolve, 5));
  settle.reject(new Error('worker failed'));
  await assert.rejects(rejected, /worker failed/);
  assert.equal(starts, 2);
  assert.deepEqual(Object.keys(page.context).filter((key) => key.startsWith('__GBDRAW_TEST_EVALUATION_')), []);
});
