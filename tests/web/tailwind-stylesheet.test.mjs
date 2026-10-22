// OV-910: the page links the precompiled Tailwind stylesheet instead of compiling
// Tailwind in the browser. The tracked css/tailwind.css must equal a fresh run of
// tools/generate_web_tailwind_css.mjs, so a class added to the page, a module, or
// the hosted-build markup, or a changed component rule, fails here until the
// stylesheet is regenerated.
import assert from 'node:assert/strict';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

import {
  classSourcePaths,
  collectMarkupClasses,
  extractClassCandidates,
  generateTailwindStylesheet,
  TAILWIND_COMPILER_PATH,
  TAILWIND_OUTPUT_PATH
} from '../../tools/generate_web_tailwind_css.mjs';

const repositoryFile = (path) => readFileSync(new URL(`../../${path}`, import.meta.url), 'utf8');

test('the tracked Tailwind stylesheet equals a fresh generator run', async () => {
  const fresh = await generateTailwindStylesheet();
  assert.ok(
    repositoryFile(TAILWIND_OUTPUT_PATH) === fresh,
    `${TAILWIND_OUTPUT_PATH} is stale; run node tools/generate_web_tailwind_css.mjs`
  );
});

test('the page links the stylesheet after its own styles and ships no Tailwind compiler', () => {
  const html = repositoryFile('gbdraw/web/index.html');
  const head = html.slice(0, html.indexOf('</head>'));
  for (const path of classSourcePaths().filter((source) => source.startsWith('gbdraw/web/'))) {
    assert.doesNotMatch(repositoryFile(path), /tailwindcss-play|text\/tailwindcss/, path);
  }
  const link = head.indexOf('<link rel="stylesheet" href="./css/tailwind.css">');
  assert.ok(link > head.lastIndexOf('</style>'), 'the stylesheet follows the last <style> of <head>');
});

// The compiler runs only here, at maintenance time: it lives under tools/, with
// its licence beside it, and not in the package the Web app ships.
test('the Tailwind compiler is a maintenance input outside the Web package', () => {
  assert.equal(TAILWIND_COMPILER_PATH, 'tools/vendor/tailwindcss/tailwindcss-play.js');
  assert.ok(existsSync(new URL(`../../${TAILWIND_COMPILER_PATH}`, import.meta.url)), TAILWIND_COMPILER_PATH);
  assert.match(repositoryFile('tools/vendor/tailwindcss/README.md'), /MIT License/);
});

test('the inventory reads every class source and keeps class tokens whole', () => {
  const sources = classSourcePaths();
  for (const path of ['gbdraw/web/index.html', 'gbdraw/web/js/app/app-setup.js', 'tools/prepare_cloudflare_pages.py', 'tools/stamp_web_build.py']) {
    assert.ok(sources.includes(path), path);
  }
  const root = mkdtempSync(join(tmpdir(), 'gbdraw-tailwind-sources-'));
  try {
    mkdirSync(join(root, 'gbdraw/web/js/app'), { recursive: true });
    for (const name of ['a.js', 'b.mjs', 'c.html', 'package.json']) writeFileSync(join(root, 'gbdraw/web/js/app', name), '');
    const modules = classSourcePaths(root).filter((path) => path.startsWith('gbdraw/web/js/'));
    assert.deepEqual(modules, ['gbdraw/web/js/app/a.js', 'gbdraw/web/js/app/b.mjs', 'gbdraw/web/js/app/c.html']);
  } finally {
    rmSync(root, { recursive: true, force: true });
  }
  const tokens = extractClassCandidates(
    `<div class="grid-cols-[minmax(0,1fr)_auto] hover:bg-blue-700" :class="{ hidden: !open, 'bg-slate-50/50': on }">`
    + "el.classList.toggle('is-active', on); const c = `px-2 ${on ? 'font-bold' : ''}`;"
  );
  for (const token of ['grid-cols-[minmax(0,1fr)_auto]', 'hover:bg-blue-700', 'hidden', 'bg-slate-50/50', 'is-active', 'px-2', 'font-bold']) {
    assert.ok(tokens.has(token), token);
  }
});

// The first batch is the markup the browser parses: a class inside a <template>,
// a <script> (x-template), or a comment arrives later, with the app.
test('the first batch holds only the classes of the parsed markup', () => {
  const html = '<body class="flex"><div @click="a > b" class="bg-white p-2"></div><!-- <i class="hidden"> -->'
    + '<template v-if="x"><b class="bg-blue-50"></b></template>'
    + '<script type="text/x-template"><span class="text-red-700"></span></script><p class="mt-1"></p></body>';
  assert.deepEqual(collectMarkupClasses(html), ['bg-white', 'flex', 'mt-1', 'p-2']);
});
