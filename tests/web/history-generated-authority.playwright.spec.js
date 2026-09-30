const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { readFileSync } = require('node:fs');
const { semantics } = require('./helpers/visual-state.cjs');
const { seeds, load, generate, switchMode, snapshot, download, popup, closeEditor } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

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
        state.adv[edit.field] = edit.value;
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
  (await import('./js/services/config.js')).buildConfigData().adv[field], field);

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
    await noDerivedStatus(page);
    await page.locator('summary[aria-label="Labels"]').click();
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

    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const restored = await inspect(page, testInfo, 'undo');
    await noDerivedStatus(page);
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
      await noDerivedStatus(fresh);
      await generate(fresh);
      expectArtifact(await inspect(fresh, testInfo, 'fresh-generate'), b);
      await noDerivedStatus(fresh);
      expect(fresh.externalRequests).toEqual([]);
    } finally { await fresh.context().close(); }

    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const redone = await inspect(page, testInfo, 'redo');
    await noDerivedStatus(page);
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


test('Live palette and its History keep the scale draft out of the committed request', async ({ browser }) => {
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
    expect(await committedScaleInterval(page)).not.toBe(12345);
    const after = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    expect(after).not.toBe(before);
    for (const action of ['undo','redo']) {
      await page.evaluate(action => window.__GBDRAW_HISTORY__[action](),action);
      await noDerivedStatus(page);
    }
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js'); state.adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    await noDerivedStatus(page);
    await page.locator('summary[aria-label="Colors"]').click();
    const paletteHelp=page.locator('[data-palette-application-help]');
    await expect(paletteHelp).toContainText('Live edit');
    await page.getByText('Palette instant preview',{exact:true}).click();
    await expect(paletteHelp).toContainText('Applies on Generate');
    const applied=await page.evaluate(()=>window.__GBDRAW_APP__.results[0].content);
    const next=await page.evaluate(()=>window.__GBDRAW_APP__.paletteNames.find(name=>name!==window.__GBDRAW_APP__.selectedPalette));
    await page.getByRole('combobox',{name:'Palette',exact:true}).selectOption(next);
    await expect(page.locator('[data-default-color-application-help]')).toContainText('Applies on Generate');
    await noDerivedStatus(page);
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
    await noDerivedStatus(page);
  } finally { await page.context().close(); }
});


test('Live global stroke keeps the scale draft out of the committed request', async ({ browser }) => {
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
    expect(await committedScaleInterval(page)).not.toBe(12345);
    expect(await page.locator('.gbdraw-preview-surface svg path[data-gbdraw-feature-id][data-gbdraw-feature-part="block"]').first().getAttribute('stroke-width')).toBe('2');
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');state.adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    await noDerivedStatus(page);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await noDerivedStatus(page);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    await noDerivedStatus(page);
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
        const value = state.adv[field];
        await window.__GBDRAW_HISTORY__.runUndoable('Scale draft', () => { state.adv[field] = draftValue; });
        return value;
      }, {field,draftValue});
      await noDerivedStatus(page);
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
        const { state } = await import('./js/state.js');state.adv[field]=prior;await window.Vue.nextTick();
      }, {field,prior});
      await noDerivedStatus(page);
      await page.evaluate(async ({field,draftValue}) => {
        const { state } = await import('./js/state.js');state.adv[field]=draftValue;await window.Vue.nextTick();
      },{field,draftValue});
      await noDerivedStatus(page);
      await expect(page.getByRole('button',{name:'Generate Diagram',exact:true})).not.toHaveAccessibleDescription(/recalculates placement/);
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
      await noDerivedStatus(fresh);
      // Record final actual rendered screens for visual inspection.
      await fresh.screenshot({path:info.outputPath(`${mode}-pending.png`),fullPage:true});
      await generate(fresh);
      await noDerivedStatus(fresh);
      const generated=await inspect(fresh,info,`${mode}-generated`);
      expect((await snapshot(fresh)).mounted).toContain('S04 retained label');
      expect((await snapshot(fresh)).mounted).toContain('#123456');
      await fresh.evaluate(() => window.__GBDRAW_HISTORY__.undo());
      expectArtifact(await inspect(fresh,info,`${mode}-undo`),loaded);
      await noDerivedStatus(fresh);
      await fresh.evaluate(() => window.__GBDRAW_HISTORY__.redo());
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
      await noDerivedStatus(page);
      const direct=await snapshot(page);
      expect(direct.mounted).toContain('S04 direct edit retained');
      const counts=await page.evaluate(()=>[window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()]);
      await page.evaluate(()=>window.failS04Reflow());
      await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.labelReflowProcessing),{timeout:180000}).toBe(false);
      await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
      await noDerivedStatus(page);
      await page.evaluate(async mode=>{
        const {state}=await import('./js/state.js');state.adv[mode==='circular'?'scale_interval':'scale_font_size']=window.s04ScalePrior;
        await window.Vue.nextTick();
      },mode);
      await noDerivedStatus(page);
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
      await noDerivedStatus(page);
      expect((await snapshot(page)).mounted).toContain('S04 retry succeeds');
      expect(page.externalRequests).toEqual([]);
    } finally {await page.context().close();}
  });
}


for (const mode of ['linear', 'circular']) {
  test(`${mode} saved live block width survives Save and Load with a separate draft`, async ({ browser }, info) => {
    test.setTimeout(240000);
    const page = await load(browser, seeds[mode]);
    let fresh;
    try {
      await generate(page);
      await noDerivedStatus(page);
      await page.locator('input[aria-label="Block Stroke Width"]').evaluate(
        element => { element.closest('details').open = true; });
      const field = page.getByLabel('Block Stroke Width', { exact: true });
      await field.fill('2');
      await field.press('Tab');
      const blocks = '.gbdraw-preview-surface path[data-gbdraw-feature-id][data-gbdraw-feature-part="block"]';
      await expect(page.locator(blocks).first()).toHaveAttribute('stroke-width', '2');
      await noDerivedStatus(page);
      const before = await snapshot(page);
      const savedPath = info.outputPath(`${mode}-live-width.json.gz`);
      const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', savedPath)));
      expect(saved.editorState.originalSvgStroke.width).toBe(2);
      expect(saved.renderRequest).toEqual(before.request);
      fresh = await load(browser, savedPath);
      await expect(fresh.locator(blocks).first()).toHaveAttribute('stroke-width', '2');
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
        const draft = config.buildConfigData();
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
          sameDraft: JSON.stringify(config.buildConfigData()) === JSON.stringify(draft),
          // N-19: the rollback restores the admitted catalog by reference.
          sameCatalog: Boolean(catalog) && window.Vue.toRaw(state.featureCatalog.value) === catalog,
          sameHistory: JSON.stringify(counts) === JSON.stringify([
            window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()
          ]), width: state.adv.block_stroke_width
        };
      });
      expect(rollback).toEqual({ status: 'error', error: { code: 'UNKNOWN', stage: 'request-validation' },
        sameRequest: true, sameResources: true, sameDraft: true, sameCatalog: true,
        sameHistory: true, width: 2 });
      expect(await draftAdv(fresh, 'scale_interval')).toBe(12345);
      expect(await committedScaleInterval(fresh)).not.toBe(12345);
      await expect(fresh.locator(blocks).first()).toHaveAttribute('stroke-width', '2');
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
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    await expect.poll(catalogState).toEqual(retained);
    expect(errors).toEqual([]);
    // N-20: the checkpoint JSON holds no copy of the catalog.
    expect(colored).toMatchObject({ checkpoint: true, checkpointHasCatalog: false });
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
    expect(await page.evaluate(async () => [
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
