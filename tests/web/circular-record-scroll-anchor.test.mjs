// When the Circular record presentation panel appears, the Settings pane goes
// back to where it was scrolled once the panel has painted
// (`createCircularRecordDisplay` in app/circular-record-controls.js). The app
// keeps that restore as a promise, `waitForCircularRecordScrollAnchor`, so the
// Gallery capture tool awaits it instead of counting the frames it takes: a
// scroll the tool makes after the wait is never undone (HmmtDNA basic 3-1).
import assert from 'node:assert/strict';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue };
const scrollingElement = { scrollTop: 40 };
globalThis.document = { activeElement: null, scrollingElement, querySelector: () => null };
const { createCircularRecordDisplay } = await import('../../gbdraw/web/js/app/circular-record-controls.js');

const tasks = () => new Promise((resolve) => setImmediate(resolve));

test('the wait for the record panel scroll anchor settles only after the pane is put back', async () => {
  /** @type {() => void} */
  let paint = () => {};
  const presented = Vue.ref(false);
  const { waitForCircularRecordScrollAnchor } = createCircularRecordDisplay(/** @type {never} */ ({
    state: { mode: Vue.ref('circular'), activeDrawing: () => ({ form: { circular_record_selector: '' } }) },
    semanticMutationAvailable: Vue.ref(true),
    circularRecordDiscoveryState: Vue.ref({ status: 'ready', inputs: [] }),
    circularSingleRecordPresentationEnabled: presented,
    sessionOperationAvailability: () => null,
    shownDrawingOfMode: () => null,
    rebindAnnotationRecords: () => {},
    revealSetting: async () => {},
    circularTrackSlotsPanelOpen: Vue.ref(false),
    circularRecordPresentationPanel: Vue.ref(null),
    afterPaint: () => new Promise((resolve) => { paint = () => resolve(undefined); })
  }));
  await waitForCircularRecordScrollAnchor();

  presented.value = true; // the read finds one record: the panel appears
  let settled = false;
  const wait = waitForCircularRecordScrollAnchor().then(() => { settled = true; });
  scrollingElement.scrollTop = 900; // a scroll made before the panel has painted
  await tasks();
  assert.equal(settled, false);

  paint();
  await wait;
  assert.equal(scrollingElement.scrollTop, 40);
});
