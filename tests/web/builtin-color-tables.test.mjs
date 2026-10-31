import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseSpecificRules, SPECIFIC_RULE_PRESETS } from '../../gbdraw/web/js/services/file-imports.js';

// The presets are the built-in color tables of gbdraw/data/color_tables/, which
// the CLI names `-t builtin:<id>`; the Web lists and applies them as Python reads them.
test('the rule presets are the built-in color tables, in Python order', () => {
  assert.deepEqual(SPECIFIC_RULE_PRESETS.map((preset) => preset.id), ['cog', 'phrog', 'plastid', 'mitochondrion']);
  for (const preset of SPECIFIC_RULE_PRESETS) {
    assert.ok(preset.label);
    assert.ok(preset.rules.length > 0);
  }
});

test('each preset holds the rows the Web reads from the same table file', () => {
  for (const preset of SPECIFIC_RULE_PRESETS) {
    const text = readFileSync(new URL(`../../gbdraw/data/color_tables/${preset.id}.tsv`, import.meta.url), 'utf8');
    assert.deepEqual(preset.rules.map((rule) => ({ ...rule })), parseSpecificRules(text).rules, preset.id);
  }
});
