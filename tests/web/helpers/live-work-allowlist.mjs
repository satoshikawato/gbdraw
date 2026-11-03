// The R3 work allowlist (gbdraw/web/CLAUDE.md, R3; PD-OI-066): one row per
// path that shows the editor intent on a Result. The work guard in
// tests/web/live-generate-parity.playwright.spec.js runs the rows in order on
// one two-Result batch page and compares what each row runs with its entry;
// the PR random walk (tests/web/live-generate-random-walk.playwright.spec.js)
// checks that it names a step kind for each row (tests/web/helpers/live-random-walk.mjs).
import { readFileSync } from 'node:fs';
import { expect, test } from '@playwright/test';
import { evaluateWithRetainedPromise } from './app-lifecycle.cjs';
import { BATCH_FIXTURE } from './audit-browser.cjs';
import { settleLive, showResult } from './live-generate-parity.cjs';
import { download, switchMode } from './mode-transition.cjs';
import {
  addColorRule, addVisibilityRule, BATCH_0004_OFF, deleteLegendRow, displayWhileSaveHeld, generate, history,
  legendRowColor, legendRowStrokeColor, loadFeatureEdits, popupEdit, renameLegendRow, renameRow, resetSettings,
  saveWhileDisplaying, switchPalette
} from './live-generate-parity-steps.cjs';

// Label Rendering = Embedded Only, generated; the first Result is displayed.
// Records the first-record features whose labels it draws.
/** @type {string[]} */
let embeddedLabelFeatures = [];
const generateEmbeddedOnly = async (page) => {
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.label_rendering = 'embedded_only'; });
  await generate(page);
  if (await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== 0) await showResult(page, 0);
  embeddedLabelFeatures = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg = app.svgContainer.querySelector('svg');
    return app.extractedFeatures.filter((item) => String(item.locus_tag || '').startsWith('TESTA_')
      && svg.querySelector(`text[data-label-feature-id="${item.svg_id}"]`)).map((item) => item.locus_tag);
  });
  expect(embeddedLabelFeatures.length, 'Embedded Only draws two labels of the first record').toBeGreaterThan(1);
};

// D-15: the popup's Apply to all on the first feature of `type`, a palette row.
const applyToAllOnPaletteRow = (page, type, color) => evaluateWithRetainedPromise(page, async ({ featureType, value }) => {
  const app = window.__GBDRAW_APP__;
  await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => item.type === featureType), null);
  await window.Vue.nextTick();
  await app.updateClickedFeatureColor(value);
  if (app.featureStyleScopeDialog.defaultColorType !== featureType) throw new Error(`${featureType} is no palette row`);
  await app.handleFeatureStyleScopeChoice('caption');
  app.clickedFeature = null;
}, { featureType: type, value: color }).then(() => settleLive(page));
// D-15: a switch to a palette not in `avoid`, answered with Keep my colors.
const switchPaletteKeeping = (page, avoid) => evaluateWithRetainedPromise(page, async (names) => {
  const app = window.__GBDRAW_APP__;
  app.selectPalette(app.paletteNames.find((name) => !names.includes(name)));
  if (!app.paletteColorsDialog.show) throw new Error('the palette switch did not ask');
  await app.handlePaletteColorsChoice('keep');
}, avoid).then(() => settleLive(page));

// The work guard (allowlist) at the app, on a two-Result batch: each edit kind
// runs exactly the compile stages of its entry (the compile's structural
// metric, live compiles only), a Result display or a History step compiles
// once, and an edit sends exactly the worker requests of its entry, each as
// often as listed (`render` is the automatic rerender). Stages are
// `COMPILE_STAGES` in app/candidate-render.js. Every path that shows edits
// on a Result has rows here: live edits, Undo, Redo, Result display, mode
// switch, Reset Settings, Session Load, and the rerender (`render`).
export const WORK_ALLOWLIST = [
  {
    // A stroke action prepares only the saved rules (OV-198), whose matches
    // are known: no request.
    kind: 'feature stroke (popup)', stages: ['strokes'], compiles: 1, requests: [],
    run: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.setClickedFeatureStrokeColorValue('#2a9d8f');
      if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice('single');
    }),
    before: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => item.locus_tag === 'TESTA_0001'), null);
      await window.Vue.nextTick();
    }).then(() => settleLive(page))
  },
  { kind: 'Legend row stroke', stages: ['strokes'], compiles: 1, requests: [], run: (page) => legendRowStrokeColor(page, 'CDS', '#e63946') },
  // From here the drawing holds a Legend row stroke. The stroke reaches the
  // features drawn in the row's color, and an edit of the fills can change
  // them or the row's key, so every display of the fills also shows the
  // strokes (OV-1985, OV-1986).
  { kind: 'Legend row color (no rule)', stages: ['legendFills'], compiles: 1, requests: [], run: (page) => legendRowColor(page, 'tRNA', '#7b2cbf') },
  // The step restores the editor's Legend rows (their colors among them): Legend rows and fills.
  { kind: 'History step (Undo of a Legend row color)', stages: ['legend', 'legendFills'], compiles: 1, requests: [], run: (page) => history(page, 'undo') },
  // The rules' matches are prepared: the palette sends no request.
  { kind: 'palette change', stages: ['fills', 'rules', 'legendFills', 'strokes'], compiles: 1, requests: [], run: (page) => switchPalette(page, 'arctic') },
  {
    // D-15: Apply to all on a palette row sets the type's default color and
    // adds no rule; the palette watcher shows it in one compile. The choice
    // first prepares the matches of the rules that could draw the row
    // (`withTargetRules`), one rule request.
    kind: 'Apply to all on a palette row', stages: ['fills', 'rules', 'legendFills', 'strokes'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => applyToAllOnPaletteRow(page, 'tRNA', '#5e60ce')
  },
  {
    // D-15: a palette switch asks while a user default color exists; Keep is
    // one History step, which the palette watcher shows in one compile.
    kind: 'palette switch, Keep my colors', stages: ['fills', 'rules', 'legendFills', 'strokes'], compiles: 1, requests: [],
    run: (page) => switchPaletteKeeping(page, ['default', 'arctic'])
  },
  {
    kind: 'Result display after a palette change and a stroke', stages: ['legend', 'fills', 'rules', 'legendFills', 'strokes'],
    compiles: 1, requests: [], run: (page) => showResult(page, 1)
  },
  {
    // One compile for the field set that changes the rules Generate reads (the
    // value); the other steps of the add leave them as they were. The rule's
    // matches are one request; it hides features, so no rerender.
    kind: 'Feature visibility rule add', stages: ['visibility'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => addVisibilityRule(page, BATCH_0004_OFF)
  },
  {
    // The rule's caption CDS gets a second color, so Python splits the CDS
    // row in two and the automatic rerender draws the Result again (OV-43);
    // its compile is Generate's plan, and the commit's own compile shows the
    // rows it adds. Requests: the commit's rule matches and the rerender. A
    // caption written in two colors stays as written (the Legend numbers the
    // names it draws), so no caption normalization request precedes the
    // matches.
    kind: 'color rule commit', stages: ['fills', 'rules', 'legendFills', 'legend', 'strokes'], compiles: 1,
    requests: ['evaluateRules', 'render'],
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0002$', color: '#2266aa', cap: 'CDS' })
  },
  {
    kind: 'Result display, nothing changed', stages: ['legend'], compiles: 1, requests: [],
    before: (page) => showResult(page, 0), run: (page) => showResult(page, 1)
  },
  {
    // U2BFIX3 review L-A: the display a Save declined shows the Legend; the
    // fills follow in one more compile once the Save ends. The color rule
    // commit's rerender left the rule matches to prepare: the display's
    // preparation is discarded while the Save runs, so the repaint prepares
    // them again.
    kind: 'Result display declined by Save, shown when Save ends', stages: ['legend', 'fills', 'rules', 'legendFills', 'strokes'],
    compiles: 2, requests: ['evaluateRules', 'evaluateRules'],
    before: async (page) => { await showResult(page, 0); await legendRowColor(page, 'tRNA', '#2a9d8f'); },
    run: (page) => displayWhileSaveHeld(page, 1)
  },
  {
    // U2BFIX3 review M-A: Save shows the fills the declined display left
    // behind before it writes the Result (the matches are prepared).
    kind: 'Save while the displayed Result lacks paint', stages: ['legend', 'fills', 'rules', 'legendFills', 'strokes'],
    compiles: 2, requests: [],
    before: async (page) => { await showResult(page, 0); await legendRowColor(page, 'tRNA', '#e9c46a'); },
    run: async (page) => {
      await Promise.all([page.waitForEvent('download'), saveWhileDisplaying(page, 1)]);
      await settleLive(page);
    }
  },
  {
    // OV-200 (U2BFIX2 review M1): only Python decides whether the label text
    // fits, so the table asks for the rerender; its one label follow places
    // the labels too, so one reflow runs (Auto Reflow on). Requests: the
    // table and the visibility rule's matches on the catalog Generate drew;
    // the rerender's rules are Python's captions already (OV-238).
    kind: 'Feature Edits TSV load (Embedded Only)', stages: ['visibility'], compiles: 1,
    requests: ['readFeatureOverrideTable', 'evaluateRules', 'render'],
    before: async (page) => {
      await generateEmbeddedOnly(page);
      await page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = true; });
    },
    run: (page) => loadFeatureEdits(page, [[embeddedLabelFeatures[0], '', 'dup edited']])
  },
  {
    // A label a table hides needs no fit decision: no rerender (Auto Reflow
    // off). The table keeps the label text it loaded before. Requests: the
    // table, and the color and visibility rule matches on the catalog the
    // rerender drew.
    kind: 'Feature Edits TSV load, label Off (Embedded Only)', stages: ['visibility'], compiles: 1,
    requests: ['readFeatureOverrideTable', 'evaluateRules', 'evaluateRules'],
    before: (page) => page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; }),
    run: (page) => loadFeatureEdits(page, [[embeddedLabelFeatures[0], '', 'dup edited'], [embeddedLabelFeatures[1], 'off', '']])
  },
  {
    // The rerender decides the fit; its rules are Python's captions already (OV-238).
    kind: 'popup label text (Embedded Only)', stages: [], compiles: 0, requests: ['render'],
    run: (page) => popupEdit(page, embeddedLabelFeatures[0], { labelText: 'dup' })
  },
  // U3a (R14-8): the Legend edits. A row's edits show live in one compile; a
  // rule change asks for the automatic rerender only when the displayed
  // Result cannot show the rows it regroups (a split, or a row another batch
  // Result draws). Captions with one color each send no caption request. An
  // edit after a rerender or a History restore first matches the rules on the
  // features of the Result it drew (one `evaluateRules`).
  {
    kind: 'Legend rename, row without features', stages: ['legend', 'legendFills', 'strokes'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => renameLegendRow(page, 'GC content', 'GC %')
  },
  {
    // The rule takes one feature of the "other proteins" row of Result 1: a split.
    kind: 'color rule commit that splits a row', stages: ['fills', 'rules', 'legendFills', 'legend', 'strokes'], compiles: 1,
    requests: ['evaluateRules', 'render'],
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^TESTA_0005$', color: '#2a9d8f', cap: 'Five' })
  },
  {
    kind: 'Legend rename, rule row one Result draws (whole row)', stages: ['fills', 'rules', 'legendFills', 'legend', 'strokes'],
    compiles: 1, requests: ['evaluateRules'], run: (page) => renameLegendRow(page, 'Five', 'Fifth')
  },
  // The rename's History step holds the Result before and after it: Undo and
  // Redo restore the Result that shows the rows, so neither compiles nor asks.
  { kind: 'History step (Undo of the rule row rename)', stages: [], compiles: 0, requests: [], run: (page) => history(page, 'undo') },
  { kind: 'History step (Redo of the rule row rename)', stages: [], compiles: 0, requests: [], run: (page) => history(page, 'redo') },
  {
    kind: 'Legend row color, rule row', stages: ['fills', 'rules', 'legendFills', 'strokes'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => legendRowColor(page, 'Fifth', '#264653')
  },
  {
    kind: 'Rules-panel color of a captioned rule', stages: ['fills', 'rules', 'legendFills', 'strokes'], compiles: 1, requests: [],
    run: (page) => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      return app.setSpecificRuleField(app.manualSpecificRules.findIndex((rule) => rule.cap === 'Fifth'), 'color', '#e76f51');
    }).then(() => settleLive(page))
  },
  {
    // Result 2 draws the row too: Python draws it again there.
    kind: 'Legend rename, rule row both Results draw', stages: ['fills', 'rules', 'legendFills', 'legend', 'strokes'], compiles: 1,
    requests: ['render'], run: (page) => renameLegendRow(page, 'CDS [#2266aa]', 'Two')
  },
  {
    // Python writes the rules that carry the row's features under the new
    // name (synthesizeColorRules), so their matches are a request. The rules
    // also match Result 2's tRNAs, so Python draws that Result again.
    kind: 'Legend rename, palette row with features', stages: ['fills', 'rules', 'legendFills', 'legend', 'strokes'], compiles: 1,
    requests: ['evaluateRules', 'evaluateRules', 'render', 'synthesizeColorRules'],
    run: (page) => renameLegendRow(page, 'tRNA', 'transfer RNA')
  },
  // A delete, its Undo and a Restore also show the strokes: a deleted row's
  // stroke leaves its features (OV-293).
  { kind: 'Legend delete', stages: ['legend', 'strokes'], compiles: 1, requests: [], run: (page) => deleteLegendRow(page, 'GC skew (+)') },
  {
    kind: 'History step (Undo of a Legend delete)', stages: ['legend', 'legendFills', 'strokes'], compiles: 1, requests: [],
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'Legend Restore', stages: ['legend', 'legendFills', 'strokes'], compiles: 1, requests: [],
    before: (page) => deleteLegendRow(page, 'GC skew (+)'),
    run: (page) => page.evaluate(() => window.__GBDRAW_APP__.restoreAllDeletedLegendEntries()).then(() => settleLive(page))
  },
  {
    kind: 'Legend sort', stages: ['legend'], compiles: 1, requests: [],
    run: (page) => page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntries('desc')).then(() => settleLive(page))
  },
  {
    // U3b review L2: a mode switch arrives at the Result its mode drew, which
    // already shows that drawing's intent (it was left unchanged), so the
    // display compiles nothing. OV-286: the arriving Result was applied with
    // another palette than the departing one (arctic, default), so the palette
    // watcher prepares the rules once more (one `evaluateRules`).
    kind: 'mode switch back to a Result', stages: [], compiles: 0, requests: ['evaluateRules'],
    before: async (page) => {
      await switchMode(page, 'linear');
      await page.evaluate(async (text) => {
        window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'batch.gb', { type: 'text/plain', lastModified: 1000 }));
        await window.Vue.nextTick();
      }, readFileSync(BATCH_FIXTURE, 'utf8'));
      await settleLive(page);
      await generate(page);
    },
    run: (page) => switchMode(page, 'circular').then(() => settleLive(page))
  },
  {
    // A Session this page saved, loaded into it: its Results already show
    // their drawing's intent, so the Load compiles nothing and asks nothing.
    kind: 'Session Load', stages: [], compiles: 0, requests: [],
    before: (page) => download(page, 'Save Session', test.info().outputPath('work-allowlist.gbdraw-session.json'), { title: '' }),
    run: async (page) => {
      const shown = () => page.evaluate(async () => {
        const { getCommittedSvgResultRuntimeIdentity } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
        const app = window.__GBDRAW_APP__;
        return app.sessionImportPending ? null : getCommittedSvgResultRuntimeIdentity(app.results[app.selectedResultIndex]);
      });
      const before = await shown();
      await page.locator('input[accept^=".json,"]').setInputFiles(test.info().outputPath('work-allowlist.gbdraw-session.json'));
      await expect.poll(async () => { const now = await shown(); return now !== null && now !== before; }, { timeout: 180_000 }).toBe(true);
      await settleLive(page);
    }
  },
  {
    // R15-3 (OV-285): a rename onto a deleted row's caption asks; its Merge is
    // the Restore (which shows the row) and then the merge's rule commit, in
    // one History checkpoint. Each shows once: two compiles. The merge's
    // rule commit prepares its candidate rules (one request) and, as a join,
    // asks Python where the joined rows go (synthesizeColorRules); the rows
    // stay on the displayed Result, so no rerender. It runs after the Load, so
    // it leaves the rule preparation of the rows above as they were (OV-286).
    kind: 'Legend rename onto a deleted row (Restore and merge)', stages: ['fills', 'legend', 'legendFills', 'rules', 'strokes'],
    compiles: 2, requests: ['evaluateRules', 'synthesizeColorRules'],
    before: (page) => deleteLegendRow(page, 'Two'),
    run: async (page) => {
      await renameRow(page, 'Fifth', 'Two');
      expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.deletedTargetKey), 'the rename asks').not.toBe('');
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.handleLegendRenameChoice('merge'));
      await settleLive(page);
    }
  },
  {
    // Reset Settings as it is today (both drawings, one History checkpoint;
    // the drawing-scoped Reset of W3 replaces it). It shows the strokes, the
    // Legend rows and colors (OV-287) and the feature visibility (OV-290) it
    // clears in one compile; the palette watcher shows the default palette's
    // fills in a compile of its own, so a Reset that changes the palette
    // compiles twice. Removing the color rules changes a Legend source, which
    // Python draws again in one rerender (OV-43, Owner 2026-10-06), and that
    // rerender applies the live setting values Reset changes (PD-OI-102,
    // PD-OI-103). One action asks for one rerender (OV-1076).
    kind: 'Reset Settings', stages: ['fills', 'legend', 'legendFills', 'strokes', 'visibility'], compiles: 2, requests: ['render'],
    run: resetSettings
  },
  {
    // A Result that was not displayed during Reset Settings, whose bytes lack the
    // feature a Generate hid: its display asks for the one rerender that draws it (OV-2293).
    kind: 'Result display after Reset Settings', stages: ['fills', 'legend', 'legendFills', 'visibility'], compiles: 1, requests: ['render'],
    before: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        app.form.labels_mode = 'out';
        app.form.multi_record_canvas = false;
        app.adv.circular_grouping_intent = 'batch';
      });
      await generate(page);
      await showResult(page, 1);
      await popupEdit(page, { biological_feature_id: 'f5aa08bf5', record_idx: 1 }, { visibility: 'off' });
      await generate(page);
      await showResult(page, 0);
      await resetSettings(page);
    },
    run: async (page) => { await showResult(page, 1); await settleLive(page); }
  }
];
