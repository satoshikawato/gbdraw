const { test } = require('@playwright/test');
const { generate, closeEditor } = require('./helpers/mode-transition.cjs');
const { load, agree, assertNonTargetsPreserved } = require('./helpers/visual-state.cjs');

test('B-02 placement History preserves every unrelated feature', async ({ browser }, info) => {
  test.setTimeout(300000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
  try {
    await generate(page);
    const before = await agree(page, info, 'before-placement');
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.openFeatureEditorFromList(app.extractedFeatures.find(feature => feature.qualifiers.gene?.[0] === 'rps12'
        && feature.qualifiers.locus_tag?.[0] === 'NitaCp049'));
    });
    const target = await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature.svg_id);
    await page.getByLabel('Feature placement', { exact: true }).selectOption('outward');
    await closeEditor(page);
    for (const action of ['Undo', 'Redo']) {
      await page.getByRole('button', { name: action, exact: true }).click();
      const observation = await agree(page, info, action);
      await assertNonTargetsPreserved(before.completed.mounted, observation.completed.mounted, [target],
        `${action}: placement must not recolor or hide unrelated features`, true);
    }
  } finally { await page.context().close(); }
});
