import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { dirname, extname, join, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const thisFile = resolve(dirname(fileURLToPath(import.meta.url)), 'app-lifecycle-contract.test.mjs');
const roots = [
  resolve(repoRoot, 'tests/web'),
  resolve(repoRoot, 'tools'),
  resolve(repoRoot, 'docs/capture'),
  resolve(repoRoot, 'gbdraw/web/gallery/tutorials'),
  resolve(repoRoot, 'gbdraw/web/js')
];
const standaloneSources = [resolve(repoRoot, 'gbdraw/web/index.html')];
const sourceExtensions = new Set(['.cjs', '.js', '.json', '.mjs', '.py']);
const workerFieldPrefix = 'diagramGeneration' + 'Worker';
const retiredWorkerFields = [`${workerFieldPrefix}Status`, `${workerFieldPrefix}Error`];
const workerReadyField = `${workerFieldPrefix}Ready`;
const workerReadyAuthority = new RegExp([
  `(?:\\.|\\?\\.)\\s*${workerReadyField}\\b`,
  `\\b${workerReadyField}\\b\\s*(?:=|:|\\?|&&|\\|\\||===?|!==?)`,
  `\\b(?:if|while)\\s*\\([^)]*\\b${workerReadyField}\\b`,
  `\\breturn\\s+${workerReadyField}\\b`,
  `\\b(?:v-if|v-show)\\s*=\\s*["'][^"']*\\b${workerReadyField}\\b`
].join('|'), 'g');

const sourceFiles = (directory) => readdirSync(directory, { withFileTypes: true })
  .flatMap((entry) => {
    const path = join(directory, entry.name);
    if (entry.isDirectory()) return sourceFiles(path);
    if (!sourceExtensions.has(extname(entry.name))) return [];
    return resolve(path) === thisFile ? [] : [path];
  });

const lineNumberAt = (source, index) => source.slice(0, index).split('\n').length;

const extractCall = (source, openParenthesis) => {
  let depth = 1;
  let index = openParenthesis + 1;
  let quote = '';
  let triple = false;
  let escaped = false;
  let lineComment = false;
  let blockComment = false;

  while (index < source.length && depth > 0) {
    const char = source[index];
    const next = source[index + 1] || '';
    if (lineComment) {
      if (char === '\n') lineComment = false;
      index += 1;
      continue;
    }
    if (blockComment) {
      if (char === '*' && next === '/') {
        blockComment = false;
        index += 2;
      } else {
        index += 1;
      }
      continue;
    }
    if (quote) {
      if (escaped) {
        escaped = false;
        index += 1;
        continue;
      }
      if (char === '\\') {
        escaped = true;
        index += 1;
        continue;
      }
      if (triple && source.slice(index, index + 3) === quote.repeat(3)) {
        quote = '';
        triple = false;
        index += 3;
        continue;
      }
      if (!triple && char === quote) quote = '';
      index += 1;
      continue;
    }
    if (char === '/' && next === '/') {
      lineComment = true;
      index += 2;
      continue;
    }
    if (char === '/' && next === '*') {
      blockComment = true;
      index += 2;
      continue;
    }
    if (['\'', '"', '`'].includes(char)) {
      quote = char;
      triple = char !== '`' && source.slice(index, index + 3) === char.repeat(3);
      index += triple ? 3 : 1;
      continue;
    }
    if (char === '(') depth += 1;
    if (char === ')') depth -= 1;
    index += 1;
  }
  return source.slice(openParenthesis + 1, Math.max(openParenthesis + 1, index - 1));
};

const violations = [];
for (const path of [...roots.flatMap(sourceFiles), ...standaloneSources]) {
  const source = readFileSync(path, 'utf8');
  const label = relative(repoRoot, path);
  for (const retiredField of retiredWorkerFields) {
    let index = source.indexOf(retiredField);
    while (index >= 0) {
      violations.push(`${label}:${lineNumberAt(source, index)} uses retired ${retiredField}`);
      index = source.indexOf(retiredField, index + retiredField.length);
    }
  }

  const retiredDefinition = new RegExp([
    `\\b(?:def|function)\\s+wait_for_${'worker'}\\s*\\(`,
    `\\b(?:const|let|var)\\s+wait_for_${'worker'}\\s*=`
  ].join('|'), 'g');
  for (const match of source.matchAll(retiredDefinition)) {
    violations.push(
      `${label}:${lineNumberAt(source, match.index)} defines retired wait_for_${'worker'}()`
    );
  }

  for (const match of source.matchAll(workerReadyAuthority)) {
    violations.push(
      `${label}:${lineNumberAt(source, match.index)} uses ${workerReadyField} as readiness authority`
    );
  }

  const lifecycleWait = /\b(?:waitForFunction|wait_for_function)\s*\(/g;
  for (const match of source.matchAll(lifecycleWait)) {
    const openParenthesis = source.indexOf('(', match.index);
    const argumentsSource = extractCall(source, openParenthesis);
    if (argumentsSource.includes(workerReadyField)) {
      violations.push(
        `${label}:${lineNumberAt(source, match.index)} waits on ${workerReadyField}`
      );
    }
  }
}

assert.deepEqual(
  violations,
  [],
  `Browser lifecycle sources contain obsolete Worker-readiness contracts:\n${violations.join('\n')}`
);

console.log('browser app lifecycle source guard passed');
