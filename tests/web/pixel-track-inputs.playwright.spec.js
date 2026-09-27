const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { load, generate, download } = require('./helpers/mode-transition.cjs');
const { reveal, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const inspect = (page, mode) => page.evaluate(async mode => {
  const { state: s } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  return {
    svg: s.results.value[s.selectedResultIndex.value].content,
    request: getCommittedCanonicalRenderRequest(),
    geometry: JSON.parse(JSON.stringify(s.trackSlotResolvedGeometry.value)),
    slots: JSON.parse(JSON.stringify(s.adv[`${mode}_track_slots`])),
    features: s.extractedFeatures.value.length,
    rules: JSON.parse(JSON.stringify(s.manualSpecificRules.map(rule => ({ ...rule, fromFile: Boolean(rule.fromFile) })))),
    annotations: JSON.parse(JSON.stringify(s.annotationSets)),
    dispatched: window.__PIXEL_REQUESTS__?.at(-1) || null,
    runs: window.__PIXEL_REQUESTS__?.length || 0,
    processing: s.processing.value,
    completedRuns: window.__MODE_EVENTS__.filter(event => event.name === 'generate.processing-cleared').length,
    error: s.errorLog.value?.summary || null
  };
}, mode);

const openStack = async (page, mode) => {
  const button = page.getByTitle('Open Custom Track Slots', { exact: true });
  if (await button.count()) await (await reveal(button)).press('Enter');
  return page.locator(`#${mode}-custom-track-slots-panel`);
};

for (const width of [1440, 390]) {
  for (const mode of ['circular', 'linear']) {
    test(`pixel fields ${mode} keyboard/error/Generate/Session/download ${width}px`, async ({ browser }, testInfo) => {
      test.setTimeout(900000);
      const seed = mode === 'circular' ? 'tobacco-chloroplast' : 'BGC0000708-BGC0000713';
      console.log(`pixel ${mode}/${width}: load`);
      const page = await load(browser, `gbdraw/web/gallery/sessions/${seed}.gbdraw-session.json`, { width, height: 1000 });
      await page.evaluate(() => {
        window.__PIXEL_REQUESTS__ = [];
        const native = Worker.prototype.postMessage;
        Worker.prototype.postMessage = function (message, ...args) {
          if (message?.type === 'run') window.__PIXEL_REQUESTS__.push(structuredClone(message.payload.request));
          return native.call(this, message, ...args);
        };
      });
      console.log(`pixel ${mode}/${width}: open stack`);
      const panel = await openStack(page, mode);
      const wrapper = panel.locator('..');
      await expect(wrapper.getByLabel('Use custom stack', { exact: true })).toBeVisible();
      await wrapper.getByLabel('Use custom stack', { exact: true }).check();
      console.log(`pixel ${mode}/${width}: stack enabled`);
      const row = page.locator(`[data-capture="${mode}-track-slot-gc_content"]`);
      const titles = mode === 'linear' ? ['Height', 'Spacing'] : ['Inner gap', 'Outer gap'];
      const fields = mode === 'linear' ? ['height', 'spacing'] : ['inner_gap_px', 'outer_gap_px'];
      const inputs = titles.map(title => row.getByTitle(title, { exact: true }));
      await expect(row).toBeVisible();
      console.log(`pixel ${mode}/${width}: baseline Generate`);
      await generate(page);
      const source = await inspect(page, mode);
      // Focus shows the muted optional-unit placeholder; blur retains the resolved Auto hint.
      await inputs[0].fill('');
      await inputs[0].focus();
      await expect(inputs[0]).toHaveAttribute('placeholder', 'px optional');
      expect(await inputs[0].evaluate(input => getComputedStyle(input, '::placeholder').color)).toBe('rgb(148, 163, 184)');
      await inputs[0].evaluate(element => element.scrollIntoView({ block: 'center' }));
      await page.screenshot({ path: testInfo.outputPath('pixel-placeholder.png') });
      let control;
      let downloadControl;
      for (const text of ['10', '10px', '10PX', '10 px', '1e1px']) {
        for (const input of inputs) { await input.fill(text); await input.press('Tab'); }
        await expect(row.getByRole('alert')).toHaveCount(0);
        await generate(page);
        const current = await inspect(page, mode);
        const slot = current.dispatched.diagramOptions.tracks[`${mode}TrackSlots`].find(slot => slot.id === 'gc_content');
        if (mode === 'linear') {
          expect(slot.height).toEqual({ value: 10, unit: 'px' });
          expect(slot.spacing).toEqual({ value: 10, unit: 'px' });
        } else { expect(slot.innerGapPx).toBe(10); expect(slot.outerGapPx).toBe(10); }
        expect(current.features).toBe(source.features);
        expect(current.rules).toEqual(source.rules);
        expect(current.annotations).toEqual(source.annotations);
        if (control) {
          expect(current.request).toEqual(control.request);
          expect(current.dispatched).toEqual(control.dispatched);
          expect(current.geometry).toEqual(control.geometry);
          expect(current.svg).toBe(control.svg);
        } else control = current;
        if (text === '10' || text === '1e1px') {
          const pending = page.waitForEvent('download');
          await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
          const file = testInfo.outputPath(`pixel-${text}.svg`);
          await (await pending).saveAs(file);
          const bytes = await fs.readFile(file);
          if (downloadControl) expect(bytes.equals(downloadControl)).toBe(true);
          else downloadControl = bytes;
        }
      }
      // Each field reports its own constraint and retains invalid draft text.
      for (let index = 0; index < inputs.length; index += 1) {
        for (const invalid of ['px', '0x10', '1e309px', '10%', '10em', '-1px', ...(mode === 'linear' && index === 0 ? ['0'] : [])]) {
          await inputs[index].fill(invalid);
          await inputs[index].press('Tab');
          const relation = mode === 'linear' && index === 0 ? 'positive' : 'nonnegative';
          await expect(row.getByRole('alert')).toContainText(`${fields[index]} must be ${relation} finite`);
          const before = await inspect(page, mode);
          expect(before.slots.find(slot => slot.id === 'gc_content')[fields[index]]).toBe(invalid);
          await page.getByRole('button', { name: 'Generate Diagram', exact: true }).press('Enter');
          await expect.poll(async () => {
            const state = await inspect(page, mode);
            return { completedRuns: state.completedRuns, processing: state.processing };
          }, { timeout: 180000 }).toEqual({ completedRuns: before.completedRuns + 1, processing: false });
          await expect.poll(async () => (await inspect(page, mode)).error).toContain(`${fields[index]} must be ${relation} finite`);
          const failed = await inspect(page, mode);
          expect(failed.svg).toBe(control.svg);
          expect(failed.request).toEqual(control.request);
          expect(failed.geometry).toEqual(control.geometry);
          expect(failed.runs).toBe(before.runs);
          expect(failed.slots.find(slot => slot.id === 'gc_content')[fields[index]]).toBe(invalid);
          if (invalid === 'px') {
            await row.getByRole('alert').evaluate(element => element.scrollIntoView({ block: 'center' }));
            await page.screenshot({ path: testInfo.outputPath(`pixel-error-${fields[index]}.png`) });
          }
        }
        await inputs[index].fill('10PX');
      }
      // Fractional and zero values keep their physical units, with height > 0.
      await inputs[0].fill('.5px');
      await inputs[1].fill('0PX');
      await generate(page);
      const fractional = await inspect(page, mode);
      const fractionSlot = fractional.dispatched.diagramOptions.tracks[`${mode}TrackSlots`].find(slot => slot.id === 'gc_content');
      expect(mode === 'linear' ? fractionSlot.height : fractionSlot.innerGapPx).toEqual(mode === 'linear' ? { value: .5, unit: 'px' } : .5);
      expect(mode === 'linear' ? fractionSlot.spacing : fractionSlot.outerGapPx).toEqual(mode === 'linear' ? { value: 0, unit: 'px' } : 0);
      for (const input of inputs) await input.fill('10PX');
      await generate(page);
      expect((await inspect(page, mode)).svg).toBe(control.svg);
      // Save an inactive row while keeping the last Result and canonical committed request.
      await row.getByTitle('Enabled', { exact: true }).uncheck();
      const sessionPath = testInfo.outputPath('pixel.gbdraw-session.json.gz');
      const savedBytes = await download(page, 'Save Session', sessionPath);
      const saved = JSON.parse(gunzipSync(savedBytes));
      const savedSlot = saved.config.adv[`${mode}_track_slots`].find(slot => slot.id === 'gc_content');
      expect(savedSlot.enabled).toBe(false);
      if (mode === 'linear') { expect(savedSlot.height).toBe('10px'); expect(savedSlot.spacing).toBe('10px'); }
      else { expect(savedSlot.inner_gap_px).toBe('10'); expect(savedSlot.outer_gap_px).toBe('10'); }
      expect(saved.renderRequest).toEqual(control.request);
      console.log(`pixel ${mode}/${width}: Session Load`);
      const restored = await load(browser, sessionPath, { width, height: 1000 });
      const restoredState = await inspect(restored, mode);
      expect(restoredState.svg).toBe(control.svg);
      expect(restoredState.slots.find(slot => slot.id === 'gc_content').enabled).toBe(false);
      expect(restoredState.slots.find(slot => slot.id === 'gc_content')[fields[0]]).toBe(mode === 'linear' ? '10px' : '10');
      expect((await getDiagramWorkerActivity(restored)).runs).toBe(0);
      await openStack(restored, mode);
      const restoredRow = restored.locator(`[data-capture="${mode}-track-slot-gc_content"]`);
      await restoredRow.getByTitle('Disabled', { exact: true }).check();
      await generate(restored);
      expect((await inspect(restored, mode)).svg).toBe(control.svg);
      await fs.writeFile(testInfo.outputPath('pixel-evidence.json'), JSON.stringify({ seed, width, mode, sourceFeatures: source.features, controlRequest: control.request, geometry: control.geometry, savedVersion: saved.version, restoredRunsBeforeGenerate: 0, externalRequests: page.externalRequests }, null, 2));
      expect(page.externalRequests).toEqual([]);
      await restored.context().close();
      await page.context().close();
    });
  }
}
