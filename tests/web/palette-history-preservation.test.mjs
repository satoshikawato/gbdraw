import assert from 'node:assert/strict';
import { test } from 'node:test';
import { fixture } from './helpers/svg-style-fixture.mjs';

const unmatched = { feat: 'CDS', qual: 'gene', val: 'absent', color: '#abcdef' };
for (const [name, colors, rules, expected] of [
  ['palette default', { default: '#d3d3d3' }, [unmatched], '#d3d3d3'],
  ['custom default', { default: '#123456' }, [unmatched], '#123456'],
  ['explicit type', { default: '#123456', unlisted_type: '#654321' }, [unmatched], '#654321'],
  ['specific rule', { default: '#123456', unlisted_type: '#654321' }, [{ feat: 'unlisted_type', qual: 'gene', val: 'sample', color: '#aabbcc' }], '#aabbcc'],
  ['hash rule precedence', { default: '#123456' }, [{ feat: 'unlisted_type', qual: 'gene', val: 'sample', color: '#aabbcc' }, { feat: 'unlisted_type', qual: 'hash', val: 'f1', color: '#778899' }], '#778899']
]) {
  test(`specific-rule replay uses ${name}`, () => {
    const { actions, attrs } = fixture(colors, rules);
    actions.applySpecificRulesToSvg();
    assert.equal(attrs.fill, expected);
  });
}
