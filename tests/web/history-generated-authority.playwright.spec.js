const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { readFileSync } = require('node:fs');
const { semantics } = require('./helpers/visual-state.cjs');
const { seeds, load, generate, snapshot, download, popup, closeEditor } = require('./helpers/mode-transition.cjs');

const status = async (page, edit = null) => {
  const current = await page.evaluate(async edit => {
  const { getGenerationApplicationStatus } = await import('./js/services/config.js');
  const calls = [], originals = [];
  const spy = (owner, key) => {
    const original = owner[key];
    originals.push(() => { owner[key] = original; });
    owner[key] = () => { calls.push(key); throw new Error(`Unexpected Status side effect: ${key}`); };
  };
  for (const key of ['Worker', 'atob', 'btoa']) spy(window, key);
  for (const key of ['arrayBuffer', 'text']) spy(Blob.prototype, key);
  spy(crypto.subtle, 'digest'); spy(SVGElement.prototype, 'cloneNode');
  try {
    if (edit) {
      const { state } = await import('./js/state.js');
      state.adv[edit.field] = edit.value;
      await window.Vue.nextTick();
    }
    const current = getGenerationApplicationStatus();
    const feedback = window.__GBDRAW_APP__.generationApplicationFeedback;
    return { ...current, feedback, calls };
  }
  finally { originals.reverse().forEach(restore => restore()); }
  }, edit);
  const label = { clean:'Applied', pending:'Pending', unknown:'Unknown', invalid:'Invalid settings', ungenerated:'Not generated' }[current.status];
  await expect(page.locator('[data-generation-application-feedback] strong')).toHaveText(label);
  await expect(page.locator('[data-generation-application-summary] strong')).toHaveText(label);
  expect(current.feedback.label).toBe(label);
  expect(current.calls).toEqual([]);
  return current;
};

const hash = value => createHash('sha256').update(value).digest('hex');
const inspect = async (page, testInfo, name) => {
  const exported = await download(page, 'SVG', testInfo.outputPath(`${name}.svg`));
  const current = await snapshot(page);
  const editable = await page.evaluate(async () =>
    (await import('./js/services/config.js')).buildConfigData());
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
    expect(await status(page)).toMatchObject({ status:'clean', calls:[] });
    await page.locator('summary[aria-label="Labels"]').click();
    const labels = page.locator('#circular-label-mode');
    await labels.focus();
    await labels.selectOption('none');
    await labels.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Change setting');
    expect(await status(page)).toMatchObject({ status:'pending', calls:[] });
    await generate(page);
    const b = await inspect(page, testInfo, 'b');
    expect(await status(page)).toMatchObject({ status:'clean', calls:[] });
    expect(b.request).not.toEqual(a.request);
    expect(b.selected).not.toBe(a.selected);

    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const restored = await inspect(page, testInfo, 'undo');
    expect(await status(page)).toMatchObject({ status:'pending', calls:[] });
    expectArtifact(restored, a);
    expect(restored.editable).toEqual(b.editable);
    await expect(labels).toHaveValue('none');
    await page.screenshot({ path: testInfo.outputPath('undo.png') });

    const path = testInfo.outputPath('undo.gbdraw-session.json.gz');
    const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', path)));
    expect.soft(saved.renderRequest).toEqual(a.request);
    expect(saved.config.form.labels_mode).toBe('none');
    const fresh = await load(browser, path);
    try {
      const loaded = await inspect(fresh, testInfo, 'fresh-load');
      expectArtifact(loaded, a);
      expect(loaded.editable).toEqual(b.editable);
      expect(await status(fresh)).toMatchObject({ status:'pending', calls:[] });
      await generate(fresh);
      expectArtifact(await inspect(fresh, testInfo, 'fresh-generate'), b);
      expect(await status(fresh)).toMatchObject({ status:'clean', calls:[] });
      expect(fresh.externalRequests).toEqual([]);
    } finally { await fresh.context().close(); }

    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const redone = await inspect(page, testInfo, 'redo');
    expect(await status(page)).toMatchObject({ status:'clean', calls:[] });
    expectArtifact(redone, b);
    expect(redone.editable).toEqual(b.editable);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await expect(labels).toHaveValue('out');
    await generate(page);
    expectArtifact(await inspect(page, testInfo, 'undo-setting-generate'), a);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});


test('Live palette and its History retain a separate scale Pending', async ({ browser }) => {
  test.setTimeout(180000);
  const page = await load(browser);
  try {
    await generate(page);
    const before = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('Change scale', () => { state.adv.scale_interval = 12345; });
      state.paletteInstantPreviewEnabled.value = true;
      await window.__GBDRAW_HISTORY__.runUndoable('Live palette', async () => {
        state.currentColors.value = { ...state.currentColors.value, CDS:'#123456' };
        await window.Vue.nextTick();
      });
    });
    const current = await status(page);
    expect(current).toMatchObject({status:'pending',calls:[]});
    expect(current.differences.some(path=>path.includes('scale.interval'))).toBe(true);
    const after = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    expect(after).not.toBe(before);
    for (const action of ['undo','redo']) {
      await page.evaluate(action => window.__GBDRAW_HISTORY__[action](),action);
      expect(await status(page)).toMatchObject({status:'pending',calls:[]});
    }
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js'); state.adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    expect(await status(page)).toMatchObject({status:'clean',calls:[]});
    await page.locator('summary[aria-label="Colors"]').click();
    const paletteHelp=page.locator('[data-palette-application-help]');
    await expect(paletteHelp).toContainText('Live edit');
    await page.getByText('Palette instant preview',{exact:true}).click();
    await expect(paletteHelp).toContainText('Applies on Generate');
    const applied=await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content);
    const next=await page.evaluate(()=>window.__GBDRAW_APP__.paletteNames.find(name=>name!==window.__GBDRAW_APP__.selectedPalette));
    await page.getByRole('combobox',{name:'Palette',exact:true}).selectOption(next);
    await expect(page.locator('[data-default-color-application-help]')).toContainText('Applies on Generate');
    expect((await status(page)).status).toBe('pending');
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content)).toBe(applied);
    await page.evaluate(async()=>{
      const {state}=await import('./js/state.js');state.currentColors.value={...state.currentColors.value,CDS:'#345678'};
      await window.Vue.nextTick();
    });
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content)).toBe(applied);
    await page.getByText('Palette instant preview',{exact:true}).click();
    await expect(paletteHelp).toContainText('Live edit');
    await expect(page.locator('[data-default-color-application-help]')).toContainText('Live edit');
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content)).not.toBe(applied);
    expect((await status(page)).status).toBe('clean');
  } finally { await page.context().close(); }
});


test('Live global stroke advances only its applied fields and preserves scale Pending', async ({ browser }) => {
  test.setTimeout(180000);
  const page = await load(browser);
  try {
    await generate(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.adv.scale_interval = 12345;
      await window.__GBDRAW_HISTORY__.runUndoable('Live stroke', async () => {
        state.adv.block_stroke_width = 2;
        await window.Vue.nextTick();
      });
    });
    const current = await status(page);
    expect(current).toMatchObject({status:'pending',calls:[]});
    expect(current.differences.every(path=>path.includes('scale.interval'))).toBe(true);
    expect(await page.locator('.gbdraw-preview-surface svg path[data-gbdraw-feature-id][data-gbdraw-feature-part="block"]').first().getAttribute('stroke-width')).toBe('2');
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');state.adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    expect(await status(page)).toMatchObject({status:'clean',calls:[]});
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    expect((await status(page)).status).not.toBe('clean');
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    expect(await status(page)).toMatchObject({status:'pending',calls:[]});
  } finally { await page.context().close(); }
});


for (const mode of ['circular', 'linear']) {
  test(`S04 ${mode} Pending and live edits continue through Save/Load, Export, Generate and History`, async ({ browser }, info) => {
    test.setTimeout(360000);
    const page = await load(browser, seeds[mode]);
    let fresh;
    try {
      // A saved Result is displayed without generation or assuming its opaque tables match.
      await status(page);
      await generate(page);
      expect((await status(page)).status).toBe('clean');
      const field = mode === 'circular' ? 'scale_interval' : 'scale_font_size';
      const draftValue = mode === 'circular' ? 12345 : 19;
      const prior = await page.evaluate(async ({field,draftValue}) => {
        const { state } = await import('./js/state.js');
        const value = state.adv[field];
        await window.__GBDRAW_HISTORY__.runUndoable('Scale draft', () => { state.adv[field] = draftValue; });
        return value;
      }, {field,draftValue});
      expect((await status(page)).status).toBe('pending');
      await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        state.autoLabelReflowEnabled.value = false;
        state.paletteInstantPreviewEnabled.value = true;
        await window.__GBDRAW_HISTORY__.runUndoable('Live color', async () => {
          state.currentColors.value = { ...state.currentColors.value, CDS:'#123456' };
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
      expect((await status(page)).status).toBe('pending');
      const history = await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]);
      await page.evaluate(() => {
        const region = document.querySelector('[data-generation-application-announcement]');
        window.s04Announcements = [];
        new MutationObserver(() => window.s04Announcements.push(region.textContent))
          .observe(region, { childList:true, subtree:true, characterData:true });
      });
      // A changing Pending value refreshes the visible observation, not the live-region text.
      for (const value of [draftValue+1, draftValue+2, draftValue]) {
        // Keep the byte/hash/Worker/clone spies installed through Vue's render,
        // so this exercises recomputation via the UI wiring, not a cached getter.
        await status(page,{field,value});
      }
      expect(await page.evaluate(() => window.s04Announcements)).toEqual([]);
      expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(history);
      await page.evaluate(async ({field,prior}) => {
        const { state } = await import('./js/state.js');state.adv[field]=prior;await window.Vue.nextTick();
      }, {field,prior});
      expect((await status(page)).status).toBe('clean');
      await page.evaluate(async ({field,draftValue}) => {
        const { state } = await import('./js/state.js');state.adv[field]=draftValue;await window.Vue.nextTick();
      },{field,draftValue});
      expect((await status(page)).status).toBe('pending');
      expect(await page.evaluate(() => window.s04Announcements.length)).toBe(2);
      await expect(page.locator('[data-generation-application-announcement]')).toHaveAttribute('aria-live','polite');
      await expect(page.getByRole('button',{name:'Generate Diagram',exact:true})).toHaveAccessibleDescription(/recalculates placement and resets zoom.*Undo restores/);
      await expect(page.getByRole('button',{name:'Save Session',exact:true})).toHaveAccessibleDescription(/current Result.*draft.*does not Generate/);
      await expect(page.getByRole('button',{name:'SVG',exact:true})).toHaveAccessibleDescription(/Export outputs the current Result.*Neither applies pending/);
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
      expect((await status(fresh)).status).toBe('pending');
      // Record final actual rendered screens for visual inspection.
      await fresh.screenshot({path:info.outputPath(`${mode}-pending.png`),fullPage:true});
      await generate(fresh);
      expect((await status(fresh)).status).toBe('clean');
      const generated=await inspect(fresh,info,`${mode}-generated`);
      expect((await snapshot(fresh)).mounted).toContain('S04 retained label');
      expect((await snapshot(fresh)).mounted).toContain('#123456');
      await fresh.evaluate(() => window.__GBDRAW_HISTORY__.undo());
      expectArtifact(await inspect(fresh,info,`${mode}-undo`),loaded);
      expect((await status(fresh)).status).toBe('pending');
      await fresh.evaluate(() => window.__GBDRAW_HISTORY__.redo());
      const redone=await inspect(fresh,info,`${mode}-redo`);
      expect(redone.selected).toBe(generated.selected);
      expect(redone.request).toEqual(generated.request);
      expect(redone.markedMounted).toBe(true);
      // Trusted History restore does not hydrate transient label-editor bindings.
      // Compare the existing comprehensive visual oracle as well as Result bytes.
      expect(await semantics(fresh,readFileSync(info.outputPath(`${mode}-redo.svg`),'utf8')))
        .toEqual(await semantics(fresh,readFileSync(info.outputPath(`${mode}-generated.svg`),'utf8')));
      expect((await status(fresh)).status).toBe('clean');
      expect(page.externalRequests).toEqual([]);expect(fresh.externalRequests).toEqual([]);
    } finally { await page.context().close();if(fresh) await fresh.context().close(); }
  });
}


for (const mode of ['circular','linear']) {
  test(`S04 ${mode} live rerender applying/error and retry preserve separate Pending`, async ({ browser }, info) => {
    test.setTimeout(240000);
    const page=await load(browser,seeds[mode]);
    try {
      await generate(page);
      await page.evaluate(async mode => {
        const { state } = await import('./js/state.js');
        window.s04ScalePrior=state.adv[mode==='circular'?'scale_interval':'scale_font_size'];
        state.adv[mode==='circular'?'scale_interval':'scale_font_size']=mode==='circular'?12345:19;
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
      expect((await status(page)).status).toBe('pending');
      const direct=await snapshot(page);
      expect(direct.mounted).toContain('S04 direct edit retained');
      const counts=await page.evaluate(()=>[window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()]);
      await page.evaluate(()=>window.failS04Reflow());
      await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.labelReflowProcessing),{timeout:180000}).toBe(false);
      await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
      const failed=await status(page);
      expect(failed.status).toBe('pending');expect(failed.unknown).toContain('$live-render');
      await page.evaluate(async mode=>{
        const {state}=await import('./js/state.js');state.adv[mode==='circular'?'scale_interval':'scale_font_size']=window.s04ScalePrior;
        await window.Vue.nextTick();
      },mode);
      expect((await status(page)).status).toBe('unknown');
      await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
      await page.evaluate(async mode=>{
        const {state}=await import('./js/state.js');state.adv[mode==='circular'?'scale_interval':'scale_font_size']=mode==='circular'?12345:19;
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
      expect((await status(page)).status).toBe('pending');
      expect((await snapshot(page)).mounted).toContain('S04 retry succeeds');
      expect(page.externalRequests).toEqual([]);
    } finally {await page.context().close();}
  });
}
