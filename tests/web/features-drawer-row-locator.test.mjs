// P17 follow-up guard: dev-staging-only functional specs find Features drawer rows by the row's
// `title`, which the template binds to formatFeatureLocation(feat). PR CI does not run those
// specs, so when #656 changed the title to the 1-based location, the stale raw `start..end`
// locators failed only on dev staging. This PR-stage guard pins both sides of that contract:
// the template binding, and that every spec builds its row locator from the same formatter.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

const REPO_ROOT = new URL('../../', import.meta.url);
const TEMPLATE = new URL('gbdraw/web/index.html', REPO_ROOT);
const SPEC_ROOT = new URL('tests/web/', REPO_ROOT);
const ROW_TITLE = ':title="formatFeatureLocation(feat)">{{ formatFeatureLocation(feat) }}</span>';
const FORMATTER_IMPORT = "import('/gbdraw/web/js/app/feature-utils.js')";
const ROW_LOCATOR = /span\[title=/;
// A row title built from raw feature coordinates, e.g. span[title="${t.start}..${t.end}"].
const RAW_COORDINATE_LOCATOR = /span\[title="\$\{[^}]*\bstart\}\.\.\$\{[^}]*\bend\}/;

const listSpecs = (dir) => readdirSync(dir, { withFileTypes: true }).flatMap((entry) => {
  const url = new URL(entry.isDirectory() ? `${entry.name}/` : entry.name, dir);
  if (entry.isDirectory()) return entry.name === 'node_modules' ? [] : listSpecs(url);
  return entry.name.endsWith('.playwright.spec.js') ? [url] : [];
});

const rowLocatorSpecs = () => listSpecs(SPEC_ROOT)
  .map((url) => ({ path: url.pathname.slice(REPO_ROOT.pathname.length), source: readFileSync(url, 'utf8') }))
  .filter(({ source }) => ROW_LOCATOR.test(source));

export const rowLocatorViolations = (source) => {
  const violations = [];
  if (RAW_COORDINATE_LOCATOR.test(source)) violations.push('builds the row title from raw start..end');
  if (!source.includes(FORMATTER_IMPORT) || !source.includes('formatFeatureLocation(')) {
    violations.push('does not derive the row title from formatFeatureLocation');
  }
  return violations;
};

test('the Features drawer row title is the shared location formatter', () => {
  const specs = rowLocatorSpecs().map(({ path }) => path);
  assert.ok(
    readFileSync(TEMPLATE, 'utf8').includes(ROW_TITLE),
    'The Features drawer row title binding in gbdraw/web/index.html changed. PR CI does not run '
      + `these specs, which find rows by that title: ${specs.join(', ')}. Update their locators, `
      + 'run them locally, and then update this guard.'
  );
});

test('specs that find Features drawer rows by title use formatFeatureLocation', () => {
  const specs = rowLocatorSpecs();
  assert.ok(specs.length > 0, 'no spec finds Features drawer rows by title; retire this guard');
  const failures = specs
    .map(({ path, source }) => ({ path, violations: rowLocatorViolations(source) }))
    .filter(({ violations }) => violations.length > 0);
  assert.deepEqual(failures, []);
});

test('the locator check rejects a raw-coordinate row title', () => {
  const raw = 'drawer.locator(`span[title="${target.start}..${target.end}"]`)';
  assert.deepEqual(rowLocatorViolations(raw), [
    'builds the row title from raw start..end',
    'does not derive the row title from formatFeatureLocation'
  ]);
  const current = [
    "const { formatFeatureLocation } = await import('/gbdraw/web/js/app/feature-utils.js');",
    'const location = formatFeatureLocation(feature);',
    'drawer.locator(`span[title="${target.location}"]`)'
  ].join('\n');
  assert.deepEqual(rowLocatorViolations(current), []);
});
