#!/usr/bin/env node
// Generate gbdraw/web/css/tailwind.css, the Web app's Tailwind stylesheet.
//
//   node tools/generate_web_tailwind_css.mjs           write the stylesheet
//   node tools/generate_web_tailwind_css.mjs --check   exit 1 when it is stale
//
// The page loads this file instead of compiling Tailwind in the browser. The
// compiler is the Tailwind Play bundle kept in tools/vendor/tailwindcss/ (a
// maintenance input: not served, deployed, or packaged), run here in a Node
// `vm` context with the few DOM objects it reads. It receives the component rules
// of gbdraw/web/css/tailwind.input.css and every class-like token of the
// sources that put classes on the page, the input it compiled in the browser
// from the live document.
//
// Order: two utilities that set the same property on one element (a static
// `bg-white` and a bound `bg-blue-50`) are decided by their order in the CSS.
// In the browser the compiler first saw the classes of the parsed markup, then
// each class the app added later, and appended the later ones after the rules
// of the same utility. The tool keeps that order with two batches: the classes
// of the static markup of index.html, then every other token. The bindings
// should not depend on that order: tests/web/tailwind-class-order.test.mjs
// fails when an element applies two such utilities at once.
// tests/web/tailwind-stylesheet.test.mjs fails when the file differs from a
// fresh run.
import { readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join, relative } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import vm from 'node:vm';

const REPOSITORY_ROOT = fileURLToPath(new URL('..', import.meta.url));
export const TAILWIND_COMPILER_PATH = 'tools/vendor/tailwindcss/tailwindcss-play.js';
export const TAILWIND_INPUT_PATH = 'gbdraw/web/css/tailwind.input.css';
export const TAILWIND_OUTPUT_PATH = 'gbdraw/web/css/tailwind.css';
// Every source that writes class names into the page: the page and its
// templates, the modules and any markup beside them, and the markup the hosted
// build inserts.
const CLASS_SOURCE_FILES = Object.freeze([
  'gbdraw/web/index.html',
  'tools/prepare_cloudflare_pages.py',
  'tools/stamp_web_build.py'
]);
const CLASS_SOURCE_DIRECTORY = 'gbdraw/web/js';

const CLASS_SOURCE_EXTENSION = /\.(?:js|mjs|html)$/;

const listClassSources = (root, directory) => readdirSync(join(root, directory), { withFileTypes: true })
  .flatMap((entry) => {
    const path = `${directory}/${entry.name}`;
    if (entry.isDirectory()) return listClassSources(root, path);
    return entry.isFile() && CLASS_SOURCE_EXTENSION.test(entry.name) ? [path] : [];
  });

export const classSourcePaths = (root = REPOSITORY_ROOT) => [
  ...CLASS_SOURCE_FILES,
  ...listClassSources(root, CLASS_SOURCE_DIRECTORY)
].sort();

// A token is a run of characters outside quotes, markup, and code delimiters;
// a bracketed arbitrary value keeps its parentheses and commas
// (grid-cols-[minmax(0,1fr)_auto]). Each token also yields its forms without
// a trailing ':' or ',' (a Vue object key such as `hidden:`) and its parts
// between parentheses and commas. Tokens that are not Tailwind classes add
// nothing to the CSS, so the inventory may be wide; it must not miss a class.
const TOKEN = /(?:\[[^\s"'`\]]*\]|[^\s"'`<>{}=;\\[\]])+/g;
const TRAILING = /[:,.]+$/;

export const extractClassCandidates = (text) => {
  const candidates = new Set();
  const add = (token) => {
    if (!token || !/[a-z]/i.test(token)) return;
    candidates.add(token);
    const trimmed = token.replace(TRAILING, '');
    if (trimmed && trimmed !== token) candidates.add(trimmed);
  };
  for (const [token] of text.matchAll(TOKEN)) {
    add(token);
    if (!token.includes('[')) token.split(/[(),]+/).forEach(add);
  }
  return candidates;
};

export const collectClassCandidates = (root = REPOSITORY_ROOT) => {
  const candidates = new Set();
  for (const path of classSourcePaths(root)) {
    for (const token of extractClassCandidates(readFileSync(join(root, path), 'utf8'))) candidates.add(token);
  }
  return [...candidates].sort();
};

// Classes of the markup the browser parses before the app starts: the `class`
// attributes of index.html outside <script>, <style>, <template>, and comments.
const MARKUP_TAG = /<!--[\s\S]*?-->|<(\/?)([a-zA-Z][\w-]*)((?:\s+[^\s=>/"']+(?:\s*=\s*(?:"[^"]*"|'[^']*'|[^\s>"']+))?)*)\s*\/?>/g;
const CLASS_ATTRIBUTE = /(?:^|\s)class\s*=\s*"([^"]*)"/;

export const collectMarkupClasses = (html) => {
  const classes = new Set();
  let rawTextEnd = null;
  let templateDepth = 0;
  for (const match of html.matchAll(MARKUP_TAG)) {
    const [, closing, rawName, attributes] = match;
    if (!rawName) continue;
    const name = rawName.toLowerCase();
    if (rawTextEnd) {
      if (closing && name === rawTextEnd) rawTextEnd = null;
      continue;
    }
    if (name === 'template') templateDepth += closing ? -1 : 1;
    if (closing || templateDepth > 0) continue;
    if (name === 'script' || name === 'style') rawTextEnd = name;
    const classAttribute = CLASS_ATTRIBUTE.exec(attributes);
    if (classAttribute) classAttribute[1].split(/\s+/).filter(Boolean).forEach((token) => classes.add(token));
  }
  return [...classes].sort();
};

// Run the Play bundle against a document holding the component rules and one
// element whose classes grow batch by batch, as the live document did.
export const compileTailwindCss = async ({ compilerSource, componentRules, batches }) => {
  let observe = null;
  const appended = [];
  const typedStyle = { nodeType: 1, tagName: 'STYLE', textContent: componentRules, getAttribute: () => 'text/tailwindcss' };
  const classedElement = { classList: [] };
  const document = {
    documentElement: {},
    body: {},
    head: { append: (element) => appended.push(element) },
    querySelectorAll: (selector) => (selector.startsWith('style') ? [typedStyle] : [classedElement]),
    createElement: () => ({ textContent: '', isConnected: true })
  };
  class MutationObserver {
    constructor(callback) { observe ??= callback; }
    observe() {}
  }
  const context = { document, MutationObserver, console: { warn() {}, log() {}, error() {} } };
  context.window = context;
  context.self = context;
  vm.createContext(context);
  vm.runInContext(compilerSource, context, { filename: TAILWIND_COMPILER_PATH });
  if (!observe) throw new Error('The Tailwind compiler did not install its document observer.');
  for (const batch of batches) {
    classedElement.classList = batch;
    await observe([]);
  }
  if (appended.length !== 1) throw new Error(`The Tailwind compiler wrote ${appended.length} stylesheets; expected 1.`);
  return appended[0].textContent;
};

const HEADER = `/*! Generated by tools/generate_web_tailwind_css.mjs from ${TAILWIND_INPUT_PATH}. Do not edit. */\n`;

export const generateTailwindStylesheet = async (root = REPOSITORY_ROOT) => {
  const markup = collectMarkupClasses(readFileSync(join(root, 'gbdraw/web/index.html'), 'utf8'));
  const css = await compileTailwindCss({
    compilerSource: readFileSync(join(root, TAILWIND_COMPILER_PATH), 'utf8'),
    componentRules: readFileSync(join(root, TAILWIND_INPUT_PATH), 'utf8'),
    batches: [markup, collectClassCandidates(root)]
  });
  return `${HEADER}${css.trim()}\n`;
};

const main = async (argv) => {
  const check = argv.includes('--check');
  const outputPath = join(REPOSITORY_ROOT, TAILWIND_OUTPUT_PATH);
  const stylesheet = await generateTailwindStylesheet();
  if (check) {
    if (readFileSync(outputPath, 'utf8') === stylesheet) return 0;
    console.error(`${TAILWIND_OUTPUT_PATH} is stale; run node tools/generate_web_tailwind_css.mjs.`);
    return 1;
  }
  writeFileSync(outputPath, stylesheet);
  console.log(`Wrote ${relative(process.cwd(), outputPath)} (${stylesheet.length} bytes).`);
  return 0;
};

if (import.meta.url === pathToFileURL(process.argv[1] ?? '').href) {
  process.exitCode = await main(process.argv.slice(2));
}
