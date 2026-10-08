import assert from 'node:assert/strict';
import { test } from 'node:test';
import { fixture } from './helpers/svg-style-fixture.mjs';
import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import { displayedFeatureAddressing } from '../../gbdraw/web/js/services/feature-override-identity.js';

// The live fill of a feature, as the compile previews the palette and the
// specific-color rules (their matches prepared by Python, R4). A type the
// palette does not list takes Python's #d3d3d3, not the palette's `default`
// key (gbdraw/features/colors.py `default_color_map.get(type, "#d3d3d3")`).
const unmatched = { feat: 'CDS', qual: 'gene', val: 'absent', color: '#abcdef' };
for (const [name, colors, rules, expected] of [
  ['palette default', { default: '#d3d3d3' }, [unmatched], '#d3d3d3'],
  ['Python\'s fallback, not the palette default', { default: '#123456' }, [unmatched], '#d3d3d3'],
  ['explicit type', { default: '#123456', unlisted_type: '#654321' }, [unmatched], '#654321'],
  ['specific rule', { default: '#123456', unlisted_type: '#654321' }, [{ feat: 'unlisted_type', qual: 'gene', val: 'sample', color: '#aabbcc' }], '#aabbcc'],
  ['hash rule precedence', { default: '#123456' }, [{ feat: 'unlisted_type', qual: 'gene', val: 'sample', color: '#aabbcc' }, { feat: 'unlisted_type', qual: 'hash', val: 'f1', color: '#778899' }], '#778899']
]) {
  test(`specific-rule replay uses ${name}`, async () => {
    const { state, ready } = fixture(colors, rules);
    await ready;
    const operations = compileDirectEditorMutationPlan({
      catalogAdmission: displayedFeatureAddressing(state.extractedFeatures.value, ['one.svg'], 0),
      manualSpecificRules: rules,
      livePreview: { domains: ['featureFills'], paletteColors: colors, drawnContext: null }
    }).operationsByResult[0];
    assert.deepEqual(operations.featureFills, [{ renderedId: 'f1', color: expected }]);
  });
}
