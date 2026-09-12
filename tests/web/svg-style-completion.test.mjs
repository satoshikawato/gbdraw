import assert from 'node:assert/strict';
import { test } from 'node:test';
import { fixture } from './helpers/svg-style-fixture.mjs';

const unmatched = { feat: 'CDS', qual: 'gene', val: 'absent', color: '#abcdef' };
test('completed style action publishes selected Result without another user action', () => {
  const { actions, state } = fixture({ default: '#d3d3d3', unlisted_type: '#abcdef' }, [unmatched]);
  actions.applySpecificRulesToSvg();
  assert.equal(JSON.parse(state.results.value[0].content).fill, '#abcdef');
});
test('depth OFF and ON publish parent visibility at the style owner', () => {
  const { actions, state } = fixture({ default: '#d3d3d3' }, []);
  for (const enabled of [false, true]) {
    state.form.show_depth = enabled;
    actions.applyTrackVisibility();
    assert.equal(JSON.parse(state.results.value[0].content).depth.display, enabled ? undefined : 'none');
  }
});
