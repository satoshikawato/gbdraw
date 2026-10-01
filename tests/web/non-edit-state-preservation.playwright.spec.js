// G-C (Web GUI audit 2026-09-30): operations that are not edits must leave
// user-owned state (config, UI, editor, feature and group state; not the
// feature catalog) unchanged. Each row names a setup and a non-edit operation.
// A row whose operation currently changes that state names its audit ID in
// knownDefect; the row runs as test.fail(true, '<ID>') until the fixing PR
// removes the ID.
const { test, expect } = require('@playwright/test');
const {
  expectNoSilentStateChange,
  generateAndWaitForResult,
  observeNonEditOperation,
  reveal
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
  await page.evaluate(async (label) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
    await app.openFeatureEditorFromList(feature, null);
    app.clickedFeature.labelText = label;
    await app.updateClickedFeatureLabelText();
    app.clickedFeature = null;
  }, text);
  await settle(page);
};

const editRecordBLabels = (page) => page.evaluate(async () => {
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
    labels: JSON.parse(JSON.stringify(app.labelTextFeatureOverrides)),
    visibility: JSON.parse(JSON.stringify(app.labelVisibilityOverrides))
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
  await hide.check();
  await settle(page);
  await hide.uncheck();
  await settle(page);
};

const generateAgain = async (page) => {
  await generateAndWaitForResult(page);
  await settle(page);
};

// Returning to Circular refreshes the record list; wait until it holds the
// records it held before, so the snapshot does not race the refresh.
const modeRoundTrip = async (page) => {
  const records = await page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length);
  await switchMode(page, 'linear');
  await switchMode(page, 'circular');
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
    knownDefect: 'FE-01',
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
    knownDefect: 'PV-07',
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
      await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
      await settle(page);
      await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
      await settle(page);
    },
    check: (observation) => expect(observation.after.history).toEqual(observation.before.history)
  },
  {
    name: 'a mode round trip leaves user-owned state unchanged',
    setup: async (page) => {
      await openWithGenBank(page, HMMT);
      await generateAndWaitForResult(page);
      await settle(page);
    },
    operation: modeRoundTrip
  },
  {
    name: 'a mode round trip keeps the Multi-Record Canvas record order',
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
    operation: modeRoundTrip
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
