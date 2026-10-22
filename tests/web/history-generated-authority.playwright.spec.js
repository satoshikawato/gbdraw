const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { readFileSync } = require('node:fs');
const { semantics } = require('./helpers/visual-state.cjs');
const { seeds, load, generate, switchMode, snapshot, download, popup, closeEditor } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise, reveal } = require('./helpers/app-lifecycle.cjs');
const { settleLive } = require('./helpers/live-generate-parity.cjs');

// PD-OI-037 revision 2: there is no derived application status. Draft edits stay
// free of Worker, byte, digest, and SVG-clone work, and never replace the
// committed canonical request, which is the applied authority.
const noDerivedStatus = async (page, edit = null) => {
  const current = await page.evaluate(async edit => {
    const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
    const calls = [], originals = [];
    const spy = (owner, key) => {
      const original = owner[key];
      originals.push(() => { owner[key] = original; });
      owner[key] = () => { calls.push(key); throw new Error(`Unexpected draft-edit side effect: ${key}`); };
    };
    for (const key of ['Worker', 'atob', 'btoa']) spy(window, key);
    for (const key of ['arrayBuffer', 'text']) spy(Blob.prototype, key);
    spy(crypto.subtle, 'digest'); spy(SVGElement.prototype, 'cloneNode');
    try {
      if (edit) {
        const { state } = await import('./js/state.js');
        state.activeDrawing().adv[edit.field] = edit.value;
        await window.Vue.nextTick();
      }
      return { committed: getCommittedCanonicalRenderRequest(), calls };
    } finally { originals.reverse().forEach(restore => restore()); }
  }, edit);
  await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);
  await expect(page.locator('[data-generation-application-summary]')).toHaveCount(0);
  expect(current.calls).toEqual([]);
  return current;
};
const committedScaleInterval = async page => (await noDerivedStatus(page))
  .committed?.diagramOptions?.configOverrides?.['objects.scale.interval'] ?? null;
const draftAdv = (page, field) => page.evaluate(async field =>
  (await import('./js/services/config.js')).buildConfigData((await import('./js/state.js')).state.activeDrawing()).adv[field], field);

const hash = value => createHash('sha256').update(value).digest('hex');
const inspect = async (page, testInfo, name) => {
  const exported = await download(page, 'SVG', testInfo.outputPath(`${name}.svg`));
  const current = await snapshot(page);
  const editable = await page.evaluate(async () =>
    (await import('./js/services/config.js')).buildConfigData((await import('./js/state.js')).state.activeDrawing()));
  const mounted = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
    return serializeCleanSvg(state.svgContainer.value.querySelector('svg'));
  });
  const result = {
    request: current.request, editable, markedMounted: current.markedMounted,
    selected: hash(current.result), mounted: hash(mounted), exported: hash(exported)
  };
  await testInfo.attach(name, { body: JSON.stringify(result), contentType: 'application/json' });
  return result;
};
const expectArtifact = (actual, expected) => {
  expect.soft(actual.selected, 'selected Result').toBe(expected.selected);
  expect.soft(actual.mounted, 'mounted SVG').toBe(expected.mounted);
  expect.soft(actual.exported, 'exported SVG').toBe(expected.exported);
  expect.soft(actual.markedMounted).toBe(true);
  expect.soft(actual.request, 'committed canonical request').toEqual(expected.request);
};

test('Undo Generate restores request A and Result A while preserving draft B through Save and fresh Load', async ({ browser }, testInfo) => {
  test.setTimeout(360000);
  const page = await load(browser);
  try {
    await generate(page);
    const a = await inspect(page, testInfo, 'a');
    await noDerivedStatus(page);
    await page.locator('summary[aria-label="Labels"]').press('Enter');
    const labels = page.locator('#circular-label-mode');
    await labels.focus();
    await labels.selectOption('none');
    await labels.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Change setting');
    await noDerivedStatus(page);
    await generate(page);
    const b = await inspect(page, testInfo, 'b');
    await noDerivedStatus(page);
    expect(b.request).not.toEqual(a.request);
    expect(b.selected).not.toBe(a.selected);

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    const restored = await inspect(page, testInfo, 'undo');
    await noDerivedStatus(page);
    expectArtifact(restored, a);
    expect(restored.editable).toEqual(b.editable);
    await expect(labels).toHaveValue('none');
    await page.screenshot({ path: testInfo.outputPath('undo.png') });

    const path = testInfo.outputPath('undo.gbdraw-session.json.gz');
    const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', path)));
    expect.soft(saved.renderRequest).toEqual(a.request);
    expect(saved.modes[saved.ui.mode].config.form.labels_mode).toBe('none');
    const fresh = await load(browser, path);
    try {
      const loaded = await inspect(fresh, testInfo, 'fresh-load');
      expectArtifact(loaded, a);
      expect(loaded.editable).toEqual(b.editable);
      await noDerivedStatus(fresh);
      await generate(fresh);
      expectArtifact(await inspect(fresh, testInfo, 'fresh-generate'), b);
      await noDerivedStatus(fresh);
      expect(fresh.externalRequests).toEqual([]);
    } finally { await fresh.context().close(); }

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
    const redone = await inspect(page, testInfo, 'redo');
    await noDerivedStatus(page);
    expectArtifact(redone, b);
    expect(redone.editable).toEqual(b.editable);
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    await expect(labels).toHaveValue('out');
    await generate(page);
    expectArtifact(await inspect(page, testInfo, 'undo-setting-generate'), a);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});


test('Live palette and its History keep the scale draft out of the committed request', async ({ browser }) => {
  test.setTimeout(180000);
  const page = await load(browser);
  try {
    await generate(page);
    const before = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('Change scale', () => { state.activeDrawing().adv.scale_interval = 12345; });
      state.paletteInstantPreviewEnabled.value = true;
      await window.__GBDRAW_HISTORY__.runUndoable('Live palette', async () => {
        state.activeDrawing().currentColors.value = { ...state.activeDrawing().currentColors.value, CDS:'#123456' };
        await window.Vue.nextTick();
      });
    });
    expect(await committedScaleInterval(page)).not.toBe(12345);
    const after = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    expect(after).not.toBe(before);
    for (const action of ['undo','redo']) {
      await page.evaluate(action => window.__GBDRAW_HISTORY__[action](),action);
      await noDerivedStatus(page);
    }
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js'); state.activeDrawing().adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    await noDerivedStatus(page);
    await page.locator('summary[aria-label="Colors"]').press('Enter');
    const paletteHelp=page.locator('[data-palette-application-help]');
    await expect(paletteHelp).toContainText('Live edit');
    await page.getByText('Palette instant preview',{exact:true}).click();
    await expect(paletteHelp).toContainText('Applies on Generate');
    const applied=await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content);
    const next=await page.evaluate(()=>window.__GBDRAW_APP__.paletteNames.find(name=>name!==window.__GBDRAW_APP__.selectedPalette));
    const paletteSelect=page.getByRole('combobox',{name:'Palette',exact:true});
    const previous=await paletteSelect.inputValue();
    const undoCount=await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount());
    await paletteSelect.selectOption(next);
    // D-15: CDS #123456 is a user default color, so the switch asks first; the
    // select keeps the palette until a choice, and the choice is one step.
    const paletteDialog=page.getByRole('dialog',{name:'Change palette'});
    await expect(paletteDialog).toContainText('Keep your 1 changed color from the "default" palette (CDS)?');
    await expect(paletteSelect).toHaveValue(previous);
    await paletteDialog.getByRole('button',{name:"Use the palette's colors"}).click();
    await expect(paletteDialog).toHaveCount(0);
    await expect(paletteSelect).toHaveValue(next);
    expect(await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoCount+1);
    await expect(page.locator('[data-default-color-application-help]')).toContainText('Applies on Generate');
    await noDerivedStatus(page);
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content)).toBe(applied);
    await page.evaluate(async()=>{
      const {state}=await import('./js/state.js');state.activeDrawing().currentColors.value={...state.activeDrawing().currentColors.value,CDS:'#345678'};
      await window.Vue.nextTick();
    });
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content)).toBe(applied);
    await page.getByText('Palette instant preview',{exact:true}).click();
    await expect(paletteHelp).toContainText('Live edit');
    await expect(page.locator('[data-default-color-application-help]')).toContainText('Live edit');
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content)).not.toBe(applied);
    await noDerivedStatus(page);
  } finally { await page.context().close(); }
});


// PD-OI-059/PD-OI-060 (D-04, D-05) and OIC-027: Applies on Generate settings
// never rewrite the current Result; only a successful Generate applies them.
const resultContent = page => page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
const resultRewrittenWithin = (page, original, timeout) => page.waitForFunction(
  original => window.__GBDRAW_APP__.results[0]?.content !== original, original, { timeout }
).then(() => true, () => false);
const expectUnchangedResult = async (page, committed, timeout = 3000) => {
  expect(await resultRewrittenWithin(page, committed, timeout), 'Result rewritten before Generate').toBe(false);
  expect(hash(await resultContent(page))).toBe(hash(committed));
};
const svgFacts = (page, content) => page.evaluate(content => {
  const root = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const text = selector => [...root.querySelectorAll(selector)].map(node => node.textContent).join(' / ');
  const definition = root.querySelector('g[data-gbdraw-role="record-definition"]');
  const blocks = [...root.querySelectorAll('[data-gbdraw-feature-part="block"]')];
  return {
    definition: definition ? text('g[data-gbdraw-role="record-definition"] text') : '',
    definitionFontSize: Number(definition?.querySelector('text')?.getAttribute('font-size')),
    plotTitle: text('#plot_title text'),
    blockStrokeWidths: [...new Set(blocks.map(node => Number(node.getAttribute('stroke-width'))))],
    blockStrokes: [...new Set(blocks.map(node => node.getAttribute('stroke')))]
  };
}, content);

test('Circular definition settings apply on Generate and keep crop length, GC% and record label', async ({ browser }, testInfo) => {
  test.setTimeout(360000);
  const page = await load(browser);
  try {
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      Object.assign(state.activeDrawing().form, {
        multi_record_canvas: false, circular_region_start: 1000, circular_region_end: 9000,
        circular_record_label: 'Custom label', circular_record_subtitle: 'Sub'
      });
      await window.Vue.nextTick();
    });
    await generate(page);
    const committed = await resultContent(page);
    const committedRequest = (await snapshot(page)).request;
    const generated = await svgFacts(page, committed);
    expect(generated.definition).toContain('Custom label');
    expect(generated.definition).toContain('8,001 bp');
    expect(generated.definition).toContain('44.41');

    const edits = [
      ['form', 'species', 'Homo sapiens'], ['form', 'strain', 'X1'], ['form', 'plot_title', 'Audit title'],
      ['adv', 'plot_title_position', 'top'], ['adv', 'plot_title_font_size', 30], ['adv', 'def_font_size', 22],
      ['adv', 'keep_full_definition_with_plot_title', true]
    ];
    for (const [owner, field, value] of edits) {
      await page.evaluate(async ({ owner, field, value }) => {
        const { state } = await import('./js/state.js');
        state.activeDrawing()[owner][field] = value;
        await window.Vue.nextTick();
      }, { owner, field, value });
      await expectUnchangedResult(page, committed);
    }
    // The visible field keeps its History step.
    const undoCount = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    const species = page.getByLabel('Species', { exact: true });
    await species.fill('Mus musculus');
    await species.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBeGreaterThan(undoCount);
    await expectUnchangedResult(page, committed);
    expect((await snapshot(page)).request).toEqual(committedRequest);
    await expect(page.locator('[data-applies-on-generate="species-strain"]')).toBeVisible();
    await page.locator('summary[aria-label="Titles and Record Labels"]').press('Enter');
    await expect(page.locator('[data-applies-on-generate="titles"]')).toBeVisible();
    const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', testInfo.outputPath('draft.gbdraw-session.json.gz'))));
    expect(hash(saved.results[0].content)).toBe(hash(committed));
    expect(saved.modes[saved.ui.mode].config.form.species).toBe('Mus musculus');
    expect(saved.modes[saved.ui.mode].config.adv.def_font_size).toBe(22);

    await generate(page);
    const applied = await svgFacts(page, await resultContent(page));
    expect(applied.definition).toContain('Custom label');
    expect(applied.definition).toContain('8,001 bp');
    expect(applied.definition).toContain('44.41');
    expect(applied.definitionFontSize).toBe(22);
    expect(applied.plotTitle).toContain('Audit title');
    await page.screenshot({ path: testInfo.outputPath('definition-applied.png') });
  } finally { await page.context().close(); }
});

test('Global stroke settings apply on Generate and leave the Result unchanged before it', async ({ browser }, testInfo) => {
  test.setTimeout(360000);
  const page = await load(browser);
  try {
    await generate(page);
    const committed = await resultContent(page);
    const committedRequest = (await snapshot(page)).request;
    expect(await committedScaleInterval(page)).not.toBe(12345);
    const width = await reveal(page.getByLabel('Block Stroke Width', { exact: true }).first());
    for (const value of ['5', '', '-1']) {
      await width.fill(value);
      await expectUnchangedResult(page, committed);
    }
    await expect(page.locator('[data-applies-on-generate="stroke"]')).toBeVisible();
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.activeDrawing().adv.scale_interval = 12345;
      await window.__GBDRAW_HISTORY__.runUndoable('Stroke draft', async () => {
        Object.assign(state.activeDrawing().adv, {
          block_stroke_color: '#ff0000', line_stroke_width: 3,
          axis_stroke_color: '#00ff00', axis_stroke_width: 4
        });
        await window.Vue.nextTick();
      });
    });
    await expectUnchangedResult(page, committed);
    expect(await committedScaleInterval(page)).not.toBe(12345);
    for (const action of ['undo', 'redo']) {
      await page.evaluate(action => window.__GBDRAW_HISTORY__[action](), action);
      await noDerivedStatus(page);
      await expectUnchangedResult(page, committed, 1000);
    }
    expect((await snapshot(page)).request).toEqual(committedRequest);

    // An invalid width is still rejected by Generate, which keeps the Result.
    const failed = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
    expect(failed.status).toBe('error');
    expect(hash(await resultContent(page))).toBe(hash(committed));

    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.activeDrawing().adv.scale_interval = null;
      state.activeDrawing().adv.block_stroke_width = 5;
      await window.Vue.nextTick();
    });
    await page.evaluate(() => { window.__GBDRAW_APP__.errorLog = null; });
    await generate(page);
    const applied = await svgFacts(page, await resultContent(page));
    expect(applied.blockStrokeWidths).toEqual([5]);
    expect(applied.blockStrokes).toEqual(['#ff0000']);
    await page.screenshot({ path: testInfo.outputPath('stroke-applied.png') });
  } finally { await page.context().close(); }
});


for (const mode of ['circular', 'linear']) {
  test(`S04 ${mode} draft and live edits continue through Save/Load, Export, Generate and History`, async ({ browser }, info) => {
    test.setTimeout(360000);
    const page = await load(browser, seeds[mode]);
    let fresh;
    try {
      // A saved Result is displayed without generation or assuming its opaque tables match.
      await noDerivedStatus(page);
      await generate(page);
      await noDerivedStatus(page);
      const field = mode === 'circular' ? 'scale_interval' : 'scale_font_size';
      const draftValue = mode === 'circular' ? 12345 : 19;
      const prior = await page.evaluate(async ({field,draftValue}) => {
        const { state } = await import('./js/state.js');
        const value = state.activeDrawing().adv[field];
        await window.__GBDRAW_HISTORY__.runUndoable('Scale draft', () => { state.activeDrawing().adv[field] = draftValue; });
        return value;
      }, {field,draftValue});
      await noDerivedStatus(page);
      await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        state.autoLabelReflowEnabled.value = false;
        state.paletteInstantPreviewEnabled.value = true;
        await window.__GBDRAW_HISTORY__.runUndoable('Live color', async () => {
          state.activeDrawing().currentColors.value = { ...state.activeDrawing().currentColors.value, CDS:'#123456' };
          await window.Vue.nextTick();
        });
      });
      await popup(page);
      await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('S04 retained label');
      await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
      await closeEditor(page);
      const edited = await inspect(page, info, `${mode}-live-pending`);
      expect(edited.editable.adv[field]).toBe(draftValue);
      expect((await snapshot(page)).mounted).toContain('S04 retained label');
      expect((await snapshot(page)).mounted).toContain('#123456');
      await noDerivedStatus(page);
      const history = await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]);
      await expect(page.locator('[data-generation-application-announcement]')).toHaveCount(0);
      // Draft edits keep the Result, History, and committed request unchanged.
      for (const value of [draftValue+1, draftValue+2, draftValue]) {
        // Keep the byte/hash/Worker/clone spies installed through Vue's render,
        // so this exercises recomputation via the UI wiring, not a cached getter.
        await noDerivedStatus(page,{field,value});
      }
      expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(history);
      await page.evaluate(async ({field,prior}) => {
        const { state } = await import('./js/state.js');state.activeDrawing().adv[field]=prior;await window.Vue.nextTick();
      }, {field,prior});
      await noDerivedStatus(page);
      await page.evaluate(async ({field,draftValue}) => {
        const { state } = await import('./js/state.js');state.activeDrawing().adv[field]=draftValue;await window.Vue.nextTick();
      },{field,draftValue});
      await noDerivedStatus(page);
      await expect(page.getByRole('button',{name:'Generate Diagram',exact:true})).not.toHaveAccessibleDescription(/recalculates placement/);
      await expect(page.getByRole('button',{name:'Save Session',exact:true})).toHaveAccessibleDescription(/current Result.*draft.*does not Generate/);
      await expect(page.getByRole('button',{name:'SVG',exact:true})).toHaveAccessibleDescription(/PNG, and PDF buttons download the current Result.*Neither applies pending/);
      const file=info.outputPath(`${mode}-pending.gbdraw-session.json.gz`);
      const saved=JSON.parse(gunzipSync(await download(page,'Save Session',file)));
      expectArtifact(await inspect(page,info,`${mode}-after-save`),edited);
      fresh=await load(browser,file);
      const loaded=await inspect(fresh,info,`${mode}-loaded`);
      // Import uses the existing mandatory SVG sanitizer; compare its exact output
      // and all visual semantics, without treating markup normalization as rerender.
      const admittedSaved=await fresh.evaluate(async content => {
        const { sanitizeSvgContent } = await import('./js/services/svg-sanitization.js');
        return sanitizeSvgContent(content);
      },saved.results[0].content);
      expectArtifact(loaded,{...edited,selected:hash(admittedSaved)});
      expect(await semantics(fresh,(await snapshot(fresh)).result))
        .toEqual(await semantics(fresh,saved.results[0].content));
      await noDerivedStatus(fresh);
      // Record final actual rendered screens for visual inspection.
      await fresh.screenshot({path:info.outputPath(`${mode}-pending.png`),fullPage:true});
      await generate(fresh);
      await noDerivedStatus(fresh);
      const generated=await inspect(fresh,info,`${mode}-generated`);
      expect((await snapshot(fresh)).mounted).toContain('S04 retained label');
      expect((await snapshot(fresh)).mounted).toContain('#123456');
      await evaluateWithRetainedPromise(fresh, () => window.__GBDRAW_HISTORY__.undo());
      expectArtifact(await inspect(fresh,info,`${mode}-undo`),loaded);
      await noDerivedStatus(fresh);
      await evaluateWithRetainedPromise(fresh, () => window.__GBDRAW_HISTORY__.redo());
      const redone=await inspect(fresh,info,`${mode}-redo`);
      expect(redone.selected).toBe(generated.selected);
      expect(redone.request).toEqual(generated.request);
      expect(redone.markedMounted).toBe(true);
      // Trusted History restore does not hydrate transient label-editor bindings.
      // Compare the existing comprehensive visual oracle as well as Result bytes.
      expect(await semantics(fresh,readFileSync(info.outputPath(`${mode}-redo.svg`),'utf8')))
        .toEqual(await semantics(fresh,readFileSync(info.outputPath(`${mode}-generated.svg`),'utf8')));
      await noDerivedStatus(fresh);
      expect(page.externalRequests).toEqual([]);expect(fresh.externalRequests).toEqual([]);
    } finally { await page.context().close();if(fresh) await fresh.context().close(); }
  });
}


for (const mode of ['circular','linear']) {
  test(`S04 ${mode} live rerender applying/error and retry keep the draft separate`, async ({ browser }, info) => {
    test.setTimeout(240000);
    const page=await load(browser,seeds[mode]);
    try {
      await generate(page);
      await page.evaluate(async mode => {
        const { state } = await import('./js/state.js');
        window.s04ScalePrior=state.activeDrawing().adv[mode==='circular'?'scale_interval':'scale_font_size'];
        state.activeDrawing().adv[mode==='circular'?'scale_interval':'scale_font_size']=mode==='circular'?12345:19;
        state.autoLabelReflowEnabled.value=true;
        window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse=()=>new Promise((resolve,reject)=>{
          window.failS04Reflow=()=>reject(new Error('S04 forced live rerender failure'));
        });
      },mode);
      const target=await popup(page);
      await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('S04 direct edit retained');
      await page.getByRole('button',{name:'Apply Label',exact:true}).click();
      await expect.poll(()=>page.evaluate(()=>Boolean(window.failS04Reflow)),{timeout:180000}).toBe(true);
      await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit applying');
      await noDerivedStatus(page);
      const direct=await snapshot(page);
      expect(direct.mounted).toContain('S04 direct edit retained');
      const counts=await page.evaluate(()=>[window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()]);
      await page.evaluate(()=>window.failS04Reflow());
      await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.labelReflowProcessing),{timeout:180000}).toBe(false);
      await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
      await noDerivedStatus(page);
      await page.evaluate(async mode=>{
        const {state}=await import('./js/state.js');state.activeDrawing().adv[mode==='circular'?'scale_interval':'scale_font_size']=window.s04ScalePrior;
        await window.Vue.nextTick();
      },mode);
      await noDerivedStatus(page);
      await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
      await page.evaluate(async mode=>{
        const {state}=await import('./js/state.js');state.activeDrawing().adv[mode==='circular'?'scale_interval':'scale_font_size']=mode==='circular'?12345:19;
        await window.Vue.nextTick();
      },mode);
      expect((await snapshot(page)).result).toBe(direct.result);
      expect(await page.evaluate(()=>[window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(counts);
      await closeEditor(page);
      await page.screenshot({path:info.outputPath(`${mode}-live-error.png`),fullPage:true});
      await page.evaluate(async target => {
        delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
        const app=window.__GBDRAW_APP__;
        app.openFeatureEditorFromList(app.extractedFeatures.find(feature=>feature.svg_id===target.featureId));
      },target);
      await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('S04 retry succeeds');
      await page.getByRole('button',{name:'Apply Label',exact:true}).click();
      await expect.poll(()=>page.evaluate(()=>({processing:window.__GBDRAW_APP__.labelReflowProcessing,error:window.__GBDRAW_APP__.labelReflowLastError})),{timeout:180000}).toEqual({processing:false,error:null});
      await expect(page.locator('[data-live-application-feedback]')).toHaveCount(0);
      await noDerivedStatus(page);
      expect((await snapshot(page)).mounted).toContain('S04 retry succeeds');
      expect(page.externalRequests).toEqual([]);
    } finally {await page.context().close();}
  });
}


for (const mode of ['linear', 'circular']) {
  // PD-OI-060 (D-05): the global block width applies on Generate; the applied
  // width then travels with the saved Result, separate from a later draft.
  test(`${mode} generated block width survives Save and Load with a separate draft`, async ({ browser }, info) => {
    test.setTimeout(300000);
    const page = await load(browser, seeds[mode]);
    let fresh;
    try {
      await generate(page);
      await noDerivedStatus(page);
      const generated = await resultContent(page);
      await page.locator('input[aria-label="Block Stroke Width"]').evaluate(
        element => { element.closest('details').open = true; });
      const field = page.getByLabel('Block Stroke Width', { exact: true });
      await field.fill('2');
      await field.press('Tab');
      await expectUnchangedResult(page, generated, 1000);
      await generate(page);
      const blocks = '.gbdraw-preview-surface path[data-gbdraw-feature-id][data-gbdraw-feature-part="block"]';
      await expect(page.locator(blocks).first()).toHaveAttribute('stroke-width', /^2(\.0)?$/);
      await noDerivedStatus(page);
      const before = await snapshot(page);
      const savedPath = info.outputPath(`${mode}-live-width.json.gz`);
      const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', savedPath)));
      expect(saved.editorState.originalSvgStroke.width).toBe(2);
      expect(saved.renderRequest).toEqual(before.request);
      fresh = await load(browser, savedPath);
      await expect(fresh.locator(blocks).first()).toHaveAttribute('stroke-width', /^2(\.0)?$/);
      await noDerivedStatus(fresh);
      expect((await snapshot(fresh)).request).toEqual(before.request);
      // Admission sanitizes persisted SVG; compare the actual sanitized Result,
      // rather than mistaking mandatory normalization for a render operation.
      expect(await fresh.evaluate(async content =>
        (await import('./js/services/svg-sanitization.js'))
          .sanitizeSvgContent(content), before.result)).toBe((await snapshot(fresh)).result);
      const activity = await fresh.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__);
      expect(activity?.constructions || 0).toBe(0);
      await noDerivedStatus(fresh, { field: 'scale_interval', value: 12345 });
      expect(await committedScaleInterval(fresh)).not.toBe(12345);
      const pendingPath = info.outputPath(`${mode}-pending-width.json.gz`);
      await download(fresh, 'Save Session', pendingPath);
      await fresh.context().close();
      fresh = await load(browser, pendingPath);
      expect(await draftAdv(fresh, 'scale_interval')).toBe(12345);
      expect(await committedScaleInterval(fresh)).not.toBe(12345);
      expect((await snapshot(fresh)).request).toEqual(before.request);
      const rollback = await evaluateWithRetainedPromise(fresh, async () => {
        const config = await import('./js/services/config.js');
        const { state } = await import('./js/state.js');
        const prior = config.canonicalRenderArtifactOwner.capture();
        const catalog = window.Vue.toRaw(state.featureCatalog.value);
        const draft = config.buildConfigData(state.activeDrawing());
        const counts = [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()];
        const file = new File([await (await fetch('/gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')).arrayBuffer()], 'replacement.json');
        const result = await config.importSession({ target: { files: [file], value: 'selected' } }, {
          beforePreviewMount() { throw new Error('Status round-trip rollback probe'); }
        });
        const current = config.canonicalRenderArtifactOwner.capture();
        return {
          status: result.status, error: { code: result.error?.code, stage: result.error?.stage },
          sameRequest: current.committedCanonicalSession === prior.committedCanonicalSession,
          sameResources: current.activeSessionResourceTable === prior.activeSessionResourceTable,
          sameDraft: JSON.stringify(config.buildConfigData(state.activeDrawing())) === JSON.stringify(draft),
          // N-19: the rollback restores the admitted catalog by reference.
          sameCatalog: Boolean(catalog) && window.Vue.toRaw(state.featureCatalog.value) === catalog,
          sameHistory: JSON.stringify(counts) === JSON.stringify([
            window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()
          ]), width: state.activeDrawing().adv.block_stroke_width
        };
      });
      expect(rollback).toEqual({ status: 'error', error: { code: 'UNKNOWN', stage: 'request-validation' },
        sameRequest: true, sameResources: true, sameDraft: true, sameCatalog: true,
        sameHistory: true, width: 2 });
      expect(await draftAdv(fresh, 'scale_interval')).toBe(12345);
      expect(await committedScaleInterval(fresh)).not.toBe(12345);
      await expect(fresh.locator(blocks).first()).toHaveAttribute('stroke-width', /^2(\.0)?$/);
      expect(page.externalRequests).toEqual([]);
      expect(fresh.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
      if (fresh) await fresh.context().close();
    }
  });
}

test('Unedited Circular Save and fresh Load keeps the committed request without Worker work', async ({ browser }, info) => {
  test.setTimeout(180000);
  const page = await load(browser);
  let fresh;
  try {
    await generate(page);
    await noDerivedStatus(page);
    const before = await snapshot(page);
    const path = info.outputPath('unedited-circular.json.gz');
    await download(page, 'Save Session', path);
    fresh = await load(browser, path);
    expect((await snapshot(fresh)).request).toEqual(before.request);
    await noDerivedStatus(fresh);
    expect(await fresh.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__?.constructions || 0)).toBe(0);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});


// SE-01, N-19, N-20 (R11): checkpoint Undo and Redo restore the admitted feature
// catalog by reference; the checkpoint JSON neither copies nor signs it.
test('Checkpoint Undo and Redo keep the admitted feature catalog through a mode round trip', async ({ browser }, info) => {
  test.setTimeout(300000);
  const page = await load(browser);
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  try {
    const catalogState = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const { isAdoptedFeatureCatalog } = await import('./js/services/feature-catalog.js');
      const catalog = window.Vue.toRaw(state.featureCatalog.value);
      window.__P09_CATALOG__ ??= catalog;
      return {
        same: catalog === window.__P09_CATALOG__,
        adopted: isAdoptedFeatureCatalog(catalog),
        features: state.extractedFeatures.value.length
      };
    });
    const history = () => page.evaluate(() => {
      const owner = window.__GBDRAW_HISTORY__;
      const checkpoint = owner.getCurrentCheckpoint();
      return {
        counts: [owner.getUndoCount(), owner.getRedoCount()],
        label: owner.undoLabel(),
        checkpointBytes: owner.getDiagnostics().checkpointEstimatedBytes,
        checkpoint: Boolean(checkpoint),
        checkpointHasCatalog: Boolean(checkpoint?.editorState && 'featureCatalog' in checkpoint.editorState)
      };
    });
    const retained = { same: true, adopted: true, features: 37 };
    expect(await catalogState()).toEqual(retained);
    const before = await history();
    const catalogBytes = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return JSON.stringify(state.featureCatalog.value).length * 2;
    });
    // A color that adds a legend entry is recorded as a checkpoint.
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
    });
    await expect.poll(async () => (await history()).counts).toEqual([before.counts[0] + 1, 0]);
    const colored = await history();
    await info.attach('n20-checkpoint-bytes', {
      body: JSON.stringify({ checkpointBytes: colored.checkpointBytes - before.checkpointBytes, catalogBytes }),
      contentType: 'application/json'
    });
    expect(await catalogState()).toEqual(retained);

    for (const direction of ['Undo', 'Redo', 'Undo']) {
      await page.getByRole('button', { name: direction, exact: true }).click();
      await expect.poll(async () => (await history()).counts)
        .toEqual(direction === 'Undo' ? before.counts.map((count, index) => count + index) : colored.counts);
      expect(await catalogState(), direction).toEqual(retained);
    }
    // A step that comes while the color's automatic rerender (OV-43) draws
    // makes it draw again (D-24 (OV-731)); its Results come with their own
    // catalog (R-5), and the mode switch waits for it (OIPC-C07), where dev's
    // switch dropped it. The catalog the settled Result holds is kept by reference
    // through the round trip, and nothing replaces it once the shown Result
    // settles again.
    await settleLive(page);
    const roundTripCatalog = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const catalog = window.Vue.toRaw(state.featureCatalog.value);
      window.__P09_ROUND_TRIP_CATALOG__ ??= catalog;
      return catalog === window.__P09_ROUND_TRIP_CATALOG__;
    });
    expect(await roundTripCatalog()).toBe(true);
    expect(await catalogState()).toMatchObject({ adopted: true, features: 37 });
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    expect(await roundTripCatalog(), 'the same catalog across the round trip').toBe(true);
    await settleLive(page);
    expect(await roundTripCatalog(), 'the same catalog once the Result settles').toBe(true);
    expect(await catalogState()).toMatchObject({ adopted: true, features: 37 });
    expect(errors).toEqual([]);
    // N-20: the checkpoint JSON holds no copy of the catalog.
    expect(colored).toMatchObject({ checkpoint: true, checkpointHasCatalog: false });
  } finally { await page.context().close(); }
});

// OV-136: a label rerender draws the committed Results again, so a loaded
// Session's Result keeps its saved name (not the `out` output prefix) in the
// catalog, the Export file name, and Save Session, also through Undo and Redo,
// which show the figures the step kept with no rerender (D-24 (OV-731)).
test('A label rerender keeps the loaded Result name through Undo and Redo (OV-136)', async ({ browser }, info) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    const kept = 'HmmtDNA_basic_circular';
    const names = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const catalog = window.Vue.toRaw(state.featureCatalog.value);
      (window.__OV136_CATALOGS__ ??= []).push(catalog);
      const seen = window.__OV136_CATALOGS__;
      return {
        rerendered: seen.length === 1 || !seen.slice(0, -1).includes(catalog),
        results: state.results.value.map(result => result.name),
        catalog: catalog.items.map(item => item.resultName),
        error: state.labelReflowLastError.value
      };
    });
    const exportName = async () => {
      const [svg] = await Promise.all([
        page.waitForEvent('download'),
        page.getByRole('button', { name: 'SVG', exact: true }).click()
      ]);
      return svg.suggestedFilename();
    };
    const expected = { rerendered: true, results: [kept], catalog: [kept], error: null };
    expect(await names()).toEqual(expected);
    expect(await exportName()).toBe(`${kept}.svg`);
    // A color that adds a Legend entry asks for the label rerender (OV-43).
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
    });
    for (const step of ['colored', 'Undo', 'Redo']) {
      if (step !== 'colored') await page.getByRole('button', { name: step, exact: true }).click();
      await settleLive(page);
      expect(await names(), step).toEqual({ ...expected, rerendered: step === 'colored' });
      expect(await exportName(), step).toBe(`${kept}.svg`);
    }
    const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', info.outputPath('ov136.gbdraw-session.json.gz'))));
    expect(saved.results.map(result => result.name)).toEqual([kept]);
    expect(saved.editorState.featureCatalog.items.map(item => item.resultName)).toEqual([kept]);
  } finally { await page.context().close(); }
});

// D-24 (OV-731): a color edit that adds a Legend row asks for the automatic
// rerender, which draws the row Python draws for the other CDS ("other
// proteins"). Once that rerender is shown, the edit's History step keeps its
// figure: Undo and Redo of the edit, and of the Generate before it, post no
// Worker run, and the Legend shows the restored rules each time.
const workerRuns = page => page.evaluate(() => (window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__?.instances || [])
  .reduce((total, instance) => total + (instance.runs?.length || 0), 0));
const shownLegend = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { shownLegendKeys } = await import('./js/services/legend-svg.js');
  return {
    keys: [...shownLegendKeys(state.svgContainer.value.querySelector('svg'))].sort(),
    rules: state.activeDrawing().manualSpecificRules.length,
    error: state.labelReflowLastError.value
  };
});
const colorFirstCds = page => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  return app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
});
const traverse = (page, steps) => page.evaluate(async steps => {
  for (const step of steps) await window.__GBDRAW_HISTORY__[step]();
}, steps);

test('Undo and Redo of a color edit show its finished rerender without another Python run (OV-731)', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await settleLive(page);
    const loaded = await shownLegend(page);
    await generate(page);
    await settleLive(page);
    const generated = await shownLegend(page);
    expect(await colorFirstCds(page)).toBe(true);
    await settleLive(page);
    const edited = await shownLegend(page);
    expect(edited.keys).toContain('other proteins');
    expect(edited.rules).toBe(generated.rules + 1);
    const runs = await workerRuns(page);
    for (const [step, expected] of [['undo', generated], ['undo', loaded], ['redo', generated], ['redo', edited]]) {
      await traverse(page, [step]);
      await settleLive(page);
      expect(await shownLegend(page), step).toEqual(expected);
      expect(await workerRuns(page), `${step} posts no Worker run`).toBe(runs);
    }
  } finally { await page.context().close(); }
});

// D-24 (OV-731): an Undo or Redo that comes while the edit's rerender is still
// drawing still draws the restored rules; no figure of the replaced rules is
// shown. The rerender the Redo asks for is then the step's figure.
test('Undo and Redo while a color edit is still being drawn draw the restored rules (OV-731)', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await settleLive(page);
    const generated = await shownLegend(page);
    const before = await workerRuns(page);
    // The edit and its Undo in one task: the Undo comes once the edit's
    // rerender has posted its Worker run and is still drawing.
    expect(await evaluateWithRetainedPromise(page, async before => {
      const app = window.__GBDRAW_APP__;
      const runs = () => (window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__?.instances || [])
        .reduce((total, instance) => total + (instance.runs?.length || 0), 0);
      const applied = await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
      while (runs() === before && app.labelReflowProcessing) await new Promise(resolve => setTimeout(resolve, 5));
      const posted = runs() - before;
      const drawing = app.labelReflowProcessing;
      await window.__GBDRAW_HISTORY__.undo();
      return { applied, posted, drawing };
    }, before)).toEqual({ applied: true, posted: 1, drawing: true });
    await settleLive(page);
    expect(await shownLegend(page), 'Undo while drawing').toEqual(generated);
    // The edit's run, and the run the Undo asks for: the rerender drawing the
    // edit is not admitted and draws the restored rules (OV-351).
    expect(await workerRuns(page), 'Undo while drawing draws again').toBe(before + 2);
    // Redo asks for the rerender again; Undo and Redo come while it draws.
    const redoing = await workerRuns(page);
    await traverse(page, ['redo', 'undo', 'redo']);
    await settleLive(page);
    const redone = await shownLegend(page);
    expect(redone.keys, 'Redo while drawing').toContain('other proteins');
    expect(redone).toMatchObject({ rules: generated.rules + 1, error: null });
    expect(await workerRuns(page), 'Redo while drawing draws again').toBeGreaterThan(redoing);
    const runs = await workerRuns(page);
    for (const [step, expected] of [['undo', generated], ['redo', redone]]) {
      await traverse(page, [step]);
      await settleLive(page);
      expect(await shownLegend(page), `${step} after the rerender`).toEqual(expected);
      expect(await workerRuns(page), `${step} after the rerender posts no Worker run`).toBe(runs);
    }
  } finally { await page.context().close(); }
});

// D-36: a rerender keeps a figure only when it drew the current paint. The
// second of two quick color edits, each adding a Legend row, waits for the
// first edit's rerender and asks for its own; History keeps that one, so Redo
// shows the row Python draws for the second edit ("other proteins").
// OV-872: the Legend rows the first edit's rerender draws while the second
// edit waits are no step of their own: each edit is one step.
test('A color edit made while another edit is drawn keeps its own rerender (OV-731)', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await settleLive(page);
    const generated = await shownLegend(page);
    const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    expect(await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const first = await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'tRNA'), '#654321');
      const drawing = app.labelReflowProcessing;
      const second = await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
      return { first, drawing, second };
    })).toEqual({ first: true, drawing: true, second: true });
    await settleLive(page);
    const edited = await shownLegend(page);
    expect(edited.keys).toEqual(expect.arrayContaining(['other proteins', 'other tRNAs']));
    expect(edited).toMatchObject({ rules: generated.rules + 2, error: null });
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()), 'one step per color edit').toBe(steps + 2);
    const runs = await workerRuns(page);
    await traverse(page, ['undo']);
    await settleLive(page);
    const undone = await shownLegend(page);
    expect(undone.keys, 'Undo of the second edit').toContain('other tRNAs');
    expect(undone.keys, 'Undo of the second edit').not.toContain('other proteins');
    expect(undone.rules).toBe(generated.rules + 1);
    await traverse(page, ['redo']);
    await settleLive(page);
    expect(await shownLegend(page), 'Redo of the second edit').toEqual(edited);
    // The first edit's rerender was shown before the second edit committed,
    // so both sides of the second edit's step are drawn.
    expect(await workerRuns(page), 'Undo and Redo post no Worker run').toBe(runs);
    await traverse(page, ['undo', 'undo']);
    await settleLive(page);
    expect(await shownLegend(page), 'two Undos undo both edits').toEqual(generated);
  } finally { await page.context().close(); }
});

// OV-873: a rerender's Results replace the shown ones under new identities,
// and the Legend move History records is keyed by Result identity. The move
// the rerender draws again is no step of its own either.
test('A color edit made while another edit is drawn after a Legend move is one step (OV-873)', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await settleLive(page);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      app.layoutRepositionMode = true;
      await window.Vue.nextTick();
      const svg = app.svgContainer.querySelector('svg');
      const legend = svg.querySelector('[data-gbdraw-composition-role="legend"]');
      const bounds = legend.getBoundingClientRect();
      const x = bounds.left + Math.min(bounds.width / 2, 8);
      const y = bounds.top + Math.min(bounds.height / 2, 8);
      const mouse = (type, dx, dy, buttons) => new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: x + dx, clientY: y + dy, buttons, view: window
      });
      const frame = () => new Promise(resolve => requestAnimationFrame(resolve));
      legend.dispatchEvent(mouse('mousedown', 0, 0, 1));
      await frame();
      svg.dispatchEvent(mouse('mousemove', 12, -9, 1));
      await frame();
      svg.dispatchEvent(mouse('mouseup', 12, -9, 0));
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Move legend');
    await page.evaluate(() => { window.__GBDRAW_APP__.layoutRepositionMode = false; });
    await settleLive(page);
    const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    expect(await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const first = await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'tRNA'), '#654321');
      const drawing = app.labelReflowProcessing;
      const second = await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
      return { first, drawing, second };
    })).toEqual({ first: true, drawing: true, second: true });
    await settleLive(page);
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()), 'one step per color edit').toBe(steps + 2);
  } finally { await page.context().close(); }
});

// D-36: a redraw that ends without a figure (a Generate during it, then
// canceled) is still owed, so the next color step does not record the figure
// it starts from as drawn: its Undo asks for the redraw, and the Legend shows
// the restored rules' row ("other proteins").
test('Undo after a lost redraw draws the restored rules (OV-731)', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await settleLive(page);
    expect(await colorFirstCds(page)).toBe(true);
    await settleLive(page);
    const colored = await shownLegend(page);
    expect(colored.keys).toContain('other proteins');
    await traverse(page, ['undo']);
    await settleLive(page);
    // The color again; a Generate starts while its rerender draws and is canceled.
    expect(await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      const hooks = window.__GBDRAW_TEST_HOOKS__;
      const prior = hooks.onSessionLifecycleEvent;
      hooks.onSessionLifecycleEvent = event => {
        prior?.(event);
        if (event.name === 'generate.processing-published') app.cancelGeneration();
      };
      try {
        const applied = await app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'CDS'), '#123456');
        const drawing = app.labelReflowProcessing;
        const outcome = await app.runAnalysis();
        return { applied, drawing, outcome: outcome?.status };
      } finally { hooks.onSessionLifecycleEvent = prior; }
    })).toEqual({ applied: true, drawing: true, outcome: 'canceled' });
    await settleLive(page);
    expect((await shownLegend(page)).keys, 'the redraw was lost').not.toContain('other proteins');
    expect(await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      return app.setFeatureColorValue(app.extractedFeatures.find(feature => feature.type === 'tRNA'), '#654321');
    })).toBe(true);
    await settleLive(page);
    const runs = await workerRuns(page);
    await traverse(page, ['undo']);
    await settleLive(page);
    expect(await shownLegend(page), 'Undo shows the restored rules').toEqual(colored);
    expect(await workerRuns(page), 'Undo asks for the redraw').toBe(runs + 1);
  } finally { await page.context().close(); }
});

// GE-06 (D-28, PD-OI-082): Undo and Redo are busy while Generate replaces the
// artifact; draft edits stay allowed (PD-OI-051).
test('Undo and Redo are busy during Generate and leave the committed request unchanged', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    const labels = page.locator('#circular-label-mode');
    await page.locator('summary[aria-label="Labels"]').click();
    for (const value of ['none', 'both']) {
      await labels.focus();
      await labels.selectOption(value);
      await labels.press('Tab');
    }
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect(labels).toHaveValue('none');
    const observe = () => page.evaluate(async () => {
      const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
      const history = window.__GBDRAW_HISTORY__;
      return {
        counts: [history.getUndoCount(), history.getRedoCount()],
        labels: window.__GBDRAW_APP__.form.labels_mode,
        request: getCommittedCanonicalRenderRequest()
      };
    });
    const before = await observe();
    expect(before.counts[1]).toBe(1);

    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => new Promise(resolve => {
        window.__P09_RELEASE_GENERATE__ = resolve;
      });
    });
    const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(() => typeof window.__P09_RELEASE_GENERATE__), { timeout: 180000 })
      .toBe('function');
    const undo = page.getByRole('button', { name: 'Undo', exact: true });
    const redo = page.getByRole('button', { name: 'Redo', exact: true });
    await expect(undo).toBeDisabled();
    await expect(redo).toBeDisabled();
    const reason = 'Generating diagram. Retry after generation finishes.';
    await expect(page.locator('[data-session-busy-reason]')).toHaveText(reason);
    await page.keyboard.press('Control+z');
    await page.keyboard.press('Control+y');
    await page.keyboard.press('Control+Shift+z');
    expect(await evaluateWithRetainedPromise(page, async () => [
      await window.__GBDRAW_HISTORY__.undo(),
      await window.__GBDRAW_HISTORY__.redo()
    ])).toEqual([{ status: 'busy', reason }, { status: 'busy', reason }]);
    expect(await observe()).toEqual(before);

    // PD-OI-051: a draft edit during Generate is still recorded.
    await page.evaluate(() => window.__GBDRAW_HISTORY__.runUndoable('Draft during Generate', () => {
      window.__GBDRAW_APP__.form.prefix = 'during-generate';
    }));
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.prefix)).toBe('during-generate');
    expect((await observe()).counts).toEqual([before.counts[0] + 1, 0]);

    await page.evaluate(() => {
      delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
      window.__P09_RELEASE_GENERATE__();
    });
    await expect.poll(() => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
    }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
    const after = await observe();
    expect(after.counts).toEqual([before.counts[0] + 2, 0]);
    expect(after.request).not.toEqual(before.request);
    await expect(undo).toBeEnabled();
    await undo.click();
    await expect.poll(async () => (await observe()).request).toEqual(before.request);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});
