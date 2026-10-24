import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'feature-selector.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-selector-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'feature-selector.js'),
  await readFile(sourcePath, 'utf8'),
  'utf8'
);

const {
  SPECIFIC_COLOR_QUALIFIER_PRESETS,
  collectFeatureQualifierNames,
  specificColorQualifierSuggestions,
  exactRegexValue
} = await import(pathToFileURL(join(tempDir, 'feature-selector.js')));

assert.equal(exactRegexValue('YP_009725295.1'), '^YP_009725295\\.1$');

{
  const suggestions = specificColorQualifierSuggestions(
    collectFeatureQualifierNames([
      { qualifiers: { custom_tag: ['one'], color: ['#ff0000'] } },
      { selector: { qualifiers: { VendorKey: ['two'] } } }
    ]),
    [{ qual: 'rule_only' }, { qual: ' custom_tag ' }]
  );
  assert.equal(SPECIFIC_COLOR_QUALIFIER_PRESETS.includes('color'), true);
  assert.deepEqual(suggestions.slice(0, SPECIFIC_COLOR_QUALIFIER_PRESETS.length), SPECIFIC_COLOR_QUALIFIER_PRESETS);
  assert.equal(suggestions.includes('custom_tag'), true);
  assert.equal(suggestions.includes('VendorKey'), true);
  assert.equal(suggestions.includes('rule_only'), true);
  assert.equal(suggestions.filter((value) => value === 'custom_tag').length, 1);
}


console.log('feature selector tests passed');
