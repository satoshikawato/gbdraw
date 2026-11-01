import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { BUILTIN_COLOR_TABLES } from '../../gbdraw/web/js/builtin-color-tables.generated.js';
import { parseSpecificRules } from '../../gbdraw/web/js/services/file-imports.js';

// The Web's rule presets (state.js `specificRulePresets`) are the built-in color
// tables of gbdraw/data/color_tables/, which the CLI names `-t builtin:<name>`.
test('the generated tables are the built-in color tables, in Python order', () => {
  assert.deepEqual(Object.keys(BUILTIN_COLOR_TABLES), ['cog', 'phrog', 'plastid', 'mitochondrion']);
  for (const table of Object.values(BUILTIN_COLOR_TABLES)) {
    assert.ok(table.title);
    assert.ok(table.rows.length > 0);
  }
});

test('each generated table holds the rows the Web reads from the same table file', () => {
  for (const [name, table] of Object.entries(BUILTIN_COLOR_TABLES)) {
    const text = readFileSync(new URL(`../../gbdraw/data/color_tables/${name}.tsv`, import.meta.url), 'utf8');
    const rows = table.rows.map(([feat, qual, val, color, cap]) => ({ feat, qual, val, color, cap }));
    assert.deepEqual(rows, parseSpecificRules(text).rules, name);
  }
});
