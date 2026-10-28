// G-C (Web GUI audit 2026-09-30): operations that are not edits must leave
// user-owned state (config, UI, editor, feature and group state; not the
// feature catalog) unchanged. Each row names a setup and a non-edit operation.
// A row whose operation currently changes that state names its audit ID in
// knownDefect; the row runs as test.fail(true, '<ID>') until the fixing PR
// removes the ID.
const { test, expect } = require('@playwright/test');
const {
  answerChoicePrompt,
  expectNoSilentStateChange,
  generateAndWaitForResult,
  observeNonEditOperation,
  reveal,
  evaluateWithRetainedPromise
} = require('./helpers/app-lifecycle.cjs');
const {
  BATCH_FIXTURE,
  HMMT,
  openBatch,
  openWithGenBank,
  selectResult,
  settle,
  switchMode
} = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const labelsOut = () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; };

const editFirstCdsLabel = async (page, text) => {
  await evaluateWithRetainedPromise(page, async (label) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
    await app.openFeatureEditorFromList(feature, null);
    app.clickedFeature.labelText = label;
    await app.updateClickedFeatureLabelText();
    app.clickedFeature = null;
  }, text);
  await settle(page);
};

const editRecordBLabels = (page) => evaluateWithRetainedPromise(page, async () => {
  const app = window.__GBDRAW_APP__;
  const trna = app.extractedFeatures.find((feature) => feature.record_id === 'TESTB' && feature.type === 'tRNA');
  const hidden = app.extractedFeatures.find((feature) => feature.locus_tag === 'TESTB_0005');
  await app.openFeatureEditorFromList(trna, null);
  app.clickedFeature.labelText = 'EDITED_B_TRNA';
  await app.updateClickedFeatureLabelText();
  await app.openFeatureEditorFromList(hidden, null);
  app.clickedFeature.labelVisibility = 'off';
  await app.updateClickedFeatureLabelText();
  app.clickedFeature = null;
  return {
    labels: Object.values(app.featureOverrides).map((row) => row.labelText).filter((value) => value !== null),
    visibility: Object.values(app.featureOverrides).map((row) => row.labelVisibility).filter((value) => value !== null)
  };
});

const depthTsv = () => {
  const lines = ['reference_name\tposition\tdepth'];
  for (let position = 1; position <= 16569; position += 100) {
    lines.push(`NC_012920.1\t${position}\t${(10 + (position % 1000) / 50).toFixed(3)}`);
  }
  return `${lines.join('\n')}\n`;
};

const openCustomStack = async (page) => {
  const button = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
  await reveal(button);
  if (await button.getAttribute('aria-expanded') !== 'true') await button.click();
  await page.getByText('Use custom stack', { exact: true }).locator('input').check();
  await page.getByTitle('Replace the current custom stack from the simple controls; this is the only action that regenerates it').click();
  await settle(page);
};

const toggleHideGcContent = async (page) => {
  const hide = await reveal(page.getByRole('checkbox', { name: 'Hide GC Content', exact: true }));
  // The custom stack draws a GC content row, so the hide asks first.
  await hide.click();
  await answerChoicePrompt(page, 'Hide GC content?', 'hide');
  await expect(hide).toBeChecked();
  await settle(page);
  await hide.uncheck();
  await settle(page);
};

const generateAgain = async (page) => {
  await generateAndWaitForResult(page);
  await settle(page);
};

const historyCounts = (page) => page.evaluate(() => [
  window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()
]);

// Returning to Circular refreshes the record list; wait until it holds the
// records it held before, so the snapshot does not race the refresh.
// R11, PD-OI-090: a drawing switch writes no setting and records no step in
// either drawing's History.
const drawingRoundTrip = async (page) => {
  const records = await page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length);
  const counts = await historyCounts(page);
  await switchMode(page, 'linear');
  await switchMode(page, 'circular');
  expect(await historyCounts(page), 'a drawing round trip records no step').toEqual(counts);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(records);
  await settle(page);
};

const ROWS = [
  {
    name: 'Result selection round trip leaves user-owned state unchanged',
    setup: openBatch,
    operation: async (page) => {
      await selectResult(page, 1);
      await selectResult(page, 0);
    }
  },
  {
    name: 'Result selection keeps label edits made on another Result',
    setup: async (page) => {
      await openBatch(page);
      await selectResult(page, 1);
      const edited = await editRecordBLabels(page);
      expect(Object.values(edited.labels)).toContain('EDITED_B_TRNA');
      expect(Object.values(edited.visibility)).toContain('off');
      await settle(page);
    },
    operation: async (page) => {
      await selectResult(page, 0);
      await selectResult(page, 1);
    }
  },
  {
    name: 'a Generate without changes leaves label edits unchanged',
    setup: async (page) => {
      await openWithGenBank(page, HMMT, labelsOut);
      await generateAndWaitForResult(page);
      await editFirstCdsLabel(page, 'NO_CHANGE_LABEL');
    },
    operation: generateAgain
  },
  {
    name: 'a Generate without changes keeps canvas padding',
    setup: async (page) => {
      await openWithGenBank(page, HMMT);
      await generateAndWaitForResult(page);
      await page.getByRole('button', { name: 'Toggle canvas padding controls' }).click();
      const right = page.locator('.preview-canvas-padding input[type=number]').nth(2);
      await right.fill('150');
      await right.press('Tab');
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.canvasPadding.right)).toBe(150);
      await settle(page);
    },
    operation: generateAgain
  },
  {
    name: 'an Undo and Redo pair restores the edited state exactly',
    setup: async (page) => {
      await openWithGenBank(page, HMMT, labelsOut);
      await generateAndWaitForResult(page);
      await editFirstCdsLabel(page, 'UNDO_REDO_LABEL');
    },
    operation: async (page) => {
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
      await settle(page);
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
      await settle(page);
    },
    check: (observation) => expect(observation.after.history).toEqual(observation.before.history)
  },
  {
    name: 'a drawing round trip that records no step leaves user-owned state unchanged',
    setup: async (page) => {
      await openWithGenBank(page, HMMT);
      await generateAndWaitForResult(page);
      await settle(page);
    },
    operation: drawingRoundTrip
  },
  {
    name: 'a drawing round trip that records no step keeps the Multi-Record Canvas record order',
    setup: async (page) => {
      await openWithGenBank(page, BATCH_FIXTURE, () => {
        Object.assign(window.__GBDRAW_APP__.form, { multi_record_canvas: true });
      });
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.adv.multi_record_positions.length)).toBe(2);
      await page.evaluate(() => window.__GBDRAW_APP__.moveCircularRecordOrderDown(0));
      await expect.poll(() => page.evaluate(() => (
        window.__GBDRAW_APP__.adv.multi_record_positions.map(({ selector }) => selector)
      ))).toEqual(['#2', '#1']);
      await generateAndWaitForResult(page);
      await settle(page);
    },
    operation: drawingRoundTrip
  },
  {
    name: 'an unrelated toggle pair leaves the custom stack unchanged',
    setup: async (page) => {
      await openWithGenBank(page, HMMT);
      await openCustomStack(page);
    },
    operation: toggleHideGcContent
  },
  {
    name: 'an unrelated toggle pair keeps a disabled Depth row disabled and in place',
    setup: async (page) => {
      await openWithGenBank(page, HMMT);
      await page.evaluate((text) => window.__GBDRAW_APP__.setCircularDepthFile(
        0, new File([text], 'sampleA.depth.tsv', { type: 'text/tab-separated-values' })
      ), depthTsv());
      await openCustomStack(page);
      const depthRow = page.getByRole('group', { name: 'Circular track slot depth', exact: true });
      await depthRow.getByTitle('Move inward').click();
      await depthRow.getByTitle('Move inward').click();
      await page.getByRole('checkbox', { name: 'Enable circular track slot depth', exact: true }).uncheck();
      await settle(page);
    },
    operation: toggleHideGcContent
  },
  {
    // G-H (R11): Undo^n then Redo^n through edits of three kinds returns to
    // the edited state and History position, and a drawing round trip keeps it.
    name: 'Undo and Redo through mixed edits, then a drawing round trip that records no step, restore the edited state',
    setup: async (page) => {
      await openWithGenBank(page, HMMT, labelsOut);
      await generateAndWaitForResult(page);
      const start = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
      await editFirstCdsLabel(page, 'UNDO_REDO_CHAIN');
      // A `ui` edit: the canvas padding.
      await page.getByRole('button', { name: 'Toggle canvas padding controls' }).click();
      const right = page.locator('.preview-canvas-padding input[type=number]').nth(2);
      await right.fill('150');
      await right.press('Tab');
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.canvasPadding.right)).toBe(150);
      const labelMode = await reveal(page.locator('#circular-label-mode'));
      await labelMode.focus();
      await page.keyboard.press('ArrowDown');
      await expect(labelMode).toHaveValue('both');
      await settle(page);
      const steps = await page.evaluate((base) => window.__GBDRAW_HISTORY__.getUndoCount() - base, start);
      expect(steps, 'one Undo step per edit').toBe(3);
    },
    operation: async (page) => {
      const [undo, redo] = await historyCounts(page);
      const padding = () => page.evaluate(() => window.__GBDRAW_APP__.canvasPadding.right);
      // The padding is the second of the three steps.
      let position = 3;
      for (const action of ['undo', 'undo', 'undo', 'redo', 'redo', 'redo']) {
        await page.evaluate((name) => window.__GBDRAW_HISTORY__[name](), action);
        await settle(page);
        position += action === 'undo' ? -1 : 1;
        const padded = (await padding()) === 150;
        expect(padded, `${action} to position ${position} restores the canvas padding`).toBe(position >= 2);
      }
      expect(await historyCounts(page), 'Redo returns to the last History position').toEqual([undo, redo]);
      await drawingRoundTrip(page);
    }
  },
  {
    // G-H (R11): the Reset Settings checkpoint restores the reset state,
    // whose track stack and multi-record layout are unset, exactly.
    name: 'an Undo and Redo pair of Reset Settings restores the reset state exactly',
    setup: async (page) => {
      await openWithGenBank(page, HMMT, labelsOut);
      await generateAndWaitForResult(page);
      await editFirstCdsLabel(page, 'RESET_LABEL');
      // Reset Settings asks first; its reset action confirms.
      await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
      await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
      await settle(page);
    },
    operation: async (page) => {
      const counts = await historyCounts(page);
      for (const action of ['undo', 'redo']) {
        await page.evaluate((name) => window.__GBDRAW_HISTORY__[name](), action);
        await settle(page);
      }
      expect(await historyCounts(page), 'Redo returns to the Reset step').toEqual(counts);
    }
  },
  {
    // R17-5 (PD-OI-090 〔Reset〕): Reset this drawing is one step of the
    // shown drawing; it keeps PNG DPI, a value one per Session, and its Undo
    // restores the drawing exactly. The Linear drawing records no step.
    name: 'an Undo of Reset this drawing restores the shown drawing; the Session values and the other drawing stay',
    setup: async (page) => {
      await openWithGenBank(page, HMMT, labelsOut);
      await generateAndWaitForResult(page);
      await editFirstCdsLabel(page, 'RESET_THIS_LABEL');
      // A Session-level edit is a project step of its own (R11).
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.runProjectStep('Change setting', () => {
        window.__GBDRAW_APP__.downloadDpi = 600;
      }));
    },
    operation: async (page) => {
      await switchMode(page, 'linear');
      const linearCounts = await historyCounts(page);
      await switchMode(page, 'circular');
      await settle(page);
      const counts = await historyCounts(page);
      await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
      await answerChoicePrompt(page, 'Reset settings?', 'reset-drawing');
      await settle(page);
      expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.undoLabel(),
        window.__GBDRAW_APP__.downloadDpi])).toEqual([counts[0] + 1, 'Reset this drawing', 600]);
      await switchMode(page, 'linear');
      expect(await historyCounts(page), 'the Linear drawing records no Reset step').toEqual(linearCounts);
      await switchMode(page, 'circular');
      await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
      await settle(page);
      expect(await historyCounts(page), 'one Undo takes the Reset step').toEqual([counts[0], 1]);
    }
  },
  {
    // PD-OI-090 〔Reset〕: Cancel changes nothing and records no step.
    name: 'Cancel in the Reset dialog changes nothing and records no step',
    setup: async (page) => {
      await openWithGenBank(page, HMMT, labelsOut);
      await generateAndWaitForResult(page);
      await editFirstCdsLabel(page, 'RESET_CANCEL_LABEL');
    },
    operation: async (page) => {
      const counts = await historyCounts(page);
      await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
      const prompt = page.locator('[data-choice-prompt]');
      await expect(prompt.getByRole('heading', { name: 'Reset settings?', exact: true })).toBeVisible();
      await prompt.getByRole('button', { name: 'Cancel', exact: true }).click();
      await expect(prompt.getByRole('heading', { name: 'Reset settings?', exact: true })).toHaveCount(0);
      await settle(page);
      expect(await historyCounts(page)).toEqual(counts);
    }
  }
];

for (const row of ROWS) {
  test(row.name, async ({ page }) => {
    if (row.knownDefect) test.fail(true, row.knownDefect);
    test.setTimeout(300_000);
    await row.setup(page);
    const observation = await observeNonEditOperation(page, () => row.operation(page));
    row.check?.(observation);
    expectNoSilentStateChange(observation, { label: row.name });
  });
}
