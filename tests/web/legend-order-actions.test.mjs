import assert from 'node:assert/strict';
import test from 'node:test';

import { createLegendSortActions } from '../../gbdraw/web/js/app/legend/sort-actions.js';

const ref = (value) => ({ value });

test('move, alphabetic sort, and default order emit one document-global order intent', async () => {
  const state = {
    legendEntries: ref([
      { caption: 'Beta' },
      { caption: 'Alpha' },
      { caption: 'Only current Result' }
    ]),
    originalLegendOrder: ref(['Alpha', 'Beta', 'Only current Result'])
  };
  const requests = [];
  const actions = createLegendSortActions({
    state,
    requestLegendOrderIntent: async (intent, label) => {
      requests.push({ intent, label });
      return true;
    }
  });

  assert.equal(await actions.moveLegendEntryUp(1), true);
  assert.deepEqual(requests.at(-1).intent, {
    kind: 'order', order: ['Alpha', 'Beta', 'Only current Result']
  });
  assert.equal(await actions.moveLegendEntryDown(0), true);
  assert.deepEqual(requests.at(-1).intent.order, ['Alpha', 'Beta', 'Only current Result']);
  assert.equal(await actions.sortLegendEntries('asc'), true);
  assert.deepEqual(requests.at(-1).intent.order, ['Alpha', 'Beta', 'Only current Result']);
  assert.equal(await actions.sortLegendEntries('desc'), true);
  assert.deepEqual(requests.at(-1).intent.order, ['Only current Result', 'Beta', 'Alpha']);
  assert.equal(await actions.sortLegendEntriesByDefault(), true);
  assert.deepEqual(requests.at(-1).intent.order, ['Alpha', 'Beta', 'Only current Result']);
  assert.equal(requests.length, 5);
  requests.forEach(({ intent }) => assert.equal(intent.kind, 'order'));
});

test('legend order actions reject unavailable and degenerate requests without mutation', () => {
  const state = {
    legendEntries: ref([{ caption: 'Only' }]),
    originalLegendOrder: ref([])
  };
  let requests = 0;
  const actions = createLegendSortActions({
    state,
    requestLegendOrderIntent: () => { requests += 1; }
  });
  assert.equal(actions.moveLegendEntryUp(0), false);
  assert.equal(actions.moveLegendEntryDown(0), false);
  assert.equal(actions.sortLegendEntries('asc'), false);
  assert.equal(actions.sortLegendEntriesByDefault(), false);
  assert.equal(requests, 0);
});
