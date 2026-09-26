const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { execFileSync } = require('node:child_process');
const { load, generate, download } = require('./helpers/mode-transition.cjs');
const { reveal, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const inspect = (page, mode) => page.evaluate(async mode => {
  const { state: s } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const { getAllFeatureLegendGroups } = await import('./js/app/legend/utils.js');
  const svg = s.results.value[s.selectedResultIndex.value].content;
  const root = new DOMParser().parseFromString(svg, 'image/svg+xml').documentElement;
  return {
    svg, request: getCommittedCanonicalRenderRequest(),
    warnings: JSON.parse(JSON.stringify(s.annotationWarnings.value)),
    annotations: JSON.parse(JSON.stringify(s.annotationSets)),
    geometry: JSON.parse(JSON.stringify(s.trackSlotResolvedGeometry.value)),
    slots: JSON.parse(JSON.stringify(s.adv[`${mode}_track_slots`])),
    rules: JSON.parse(JSON.stringify(s.manualSpecificRules.map(rule => ({ ...rule, fromFile: Boolean(rule.fromFile) })))),
    legends: getAllFeatureLegendGroups(root).map(group => [...group.querySelectorAll('g[data-legend-key]')].map(entry => ({
      caption: entry.getAttribute('data-legend-key'), color: entry.querySelector('path[fill]')?.getAttribute('fill')
    }))),
    features: s.extractedFeatures.value.length,
    requestSent: window.__INTEGRATION_REQUESTS__?.at(-1),
    runs: window.__INTEGRATION_REQUESTS__?.length || 0,
    completed: window.__MODE_EVENTS__.filter(event => event.name === 'generate.processing-cleared').length,
    undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
    processing: s.processing.value, error: s.errorLog.value?.summary || null
  };
}, mode);
const drawingSignature = (page, svg) => page.evaluate(async svg => {
  // History admission sanitizes the artifact; compare every resulting drawing node
  // and attribute after that same sanitizer, including all geometry and metadata.
  const { sanitizeSvgContent } = await import('./js/services/svg-sanitization.js');
  const root = new DOMParser().parseFromString(sanitizeSvgContent(svg), 'image/svg+xml').documentElement;
  const nodes = [root, ...root.querySelectorAll('*')].map(node => [node.localName,
    [...node.attributes].map(attribute => [attribute.name, attribute.value]).sort((a, b) => a[0].localeCompare(b[0])),
    [...node.childNodes].filter(child => child.nodeType === 3).map(child => child.textContent.trim()).filter(Boolean)]);
  const digest = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(JSON.stringify(nodes)));
  return Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, '0')).join('');
}, svg);
const artifact = value => ({ svg: value.svg, warnings: value.warnings, request: value.request, geometry: value.geometry });
const history = async (page, name) => {
  await page.getByRole('button', { name, exact: true }).click();
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.restoring.value && !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
};
const colorFile = page => page.evaluate(async () => {
  const file = window.__GBDRAW_APP__.files.t_color;
  const { readFileText } = await import('./js/services/file-content-cache.js');
  return { name: file?.name || null, text: file ? await readFileText(file) : null };
});
const openStack = async (page, mode) => {
  const button = page.getByTitle('Open Custom Track Slots', { exact: true });
  if (await button.count()) await (await reveal(button)).press('Enter');
  return page.locator(`#${mode}-custom-track-slots-panel`);
};

for (const width of [1440, 390]) {
  for (const mode of ['circular', 'linear']) {
    test(`four annotation/style outcomes combine through History/Session/export ${mode} ${width}px`, async ({ browser }, testInfo) => {
      test.setTimeout(900000);
      const seed = mode === 'circular' ? 'tobacco-chloroplast' : 'BGC0000708-BGC0000713';
      const page = await load(browser, `gbdraw/web/gallery/sessions/${seed}.gbdraw-session.json`, { width, height: 1000 });
      let restored;
      const consoleText = [];
      page.on('console', message => consoleText.push(message.text()));
      try {
        await page.evaluate(() => {
          window.__INTEGRATION_REQUESTS__ = [];
          const native = Worker.prototype.postMessage;
          Worker.prototype.postMessage = function (message, ...args) {
            if (message?.type === 'run') window.__INTEGRATION_REQUESTS__.push(structuredClone(message.payload.request));
            return native.call(this, message, ...args);
          };
        });
        console.log(`integration ${mode}/${width}: loaded`);
        await generate(page);
        console.log(`integration ${mode}/${width}: baseline generated`);
        const source = await inspect(page, mode);
        const input = await page.evaluate(async () => {
          const a = window.__GBDRAW_APP__;
          const { encodeAnnotationTable, parseAnnotationTable } = await import('./js/app/annotations/table-codec.js');
          const { featureTarget } = await import('./js/app/annotations/target-actions.js');
          const { getFeatureGenerationHash } = await import('./js/app/feature-utils.js');
          const ids = a.extractedFeatures.filter(feature => feature.type === 'CDS').slice(0, 2).map(getFeatureGenerationHash);
          const sets = JSON.parse(JSON.stringify(a.annotationSets));
          if (!sets.length) sets.push(...parseAnnotationTable('set_id\tid\tmark\trecord\tstart\tend\tlabel\ncomparison_region\tanchor\tbracket\t#1\t1\t1000\tCluster start\n'));
          const record = sets[0].annotations[0].target.record;
          sets[0].annotations.push({ ...sets[0].annotations[0], id: 'integration-miss', label: 'SKIPPED-INTEGRATION-MARK', legendLabel: 'SKIPPED-INTEGRATION-LEGEND',
            target: { ...featureTarget({ selector: `hash=${ids[0]};gene=PRIVATE-INTEGRATION-MISSING` }), record } });
          const known = encodeAnnotationTable(sets);
          const auxiliary = known.trimEnd().split('\n').map((line, index) => `${line}\t${index === 0 ? 'notes\tgene_desc\tpmid' : 'PRIVATE-AUX-CELL\tPRIVATE-DESC\t12345'}`).join('\n') + '\n';
          const originalRules = a.manualSpecificRules.map(rule => [rule.feat, rule.qual, rule.val, rule.color, rule.cap].join('\t')).join('\n');
          const colors = originalRules + '\n' + ids.map((id, index) => ['CDS', 'hash', id, ['#112233', '#445566'][index], 'Integration group'].join('\t')).join('\n') + '\n';
          return { ids, auxiliary, known, colors, setId: sets[0].id };
        });
        const annotationPanel = page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
        if (await annotationPanel.getAttribute('open') === null) await annotationPanel.locator('summary').press('Enter');
        const picker = page.waitForEvent('filechooser');
        await annotationPanel.getByRole('button', { name: 'Import TSV', exact: true }).press('Enter');
        await (await picker).setFiles({ name: 'integration-annotations.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(input.auxiliary) });
        await expect.poll(() => annotationPanel.locator('input[type=file]').inputValue()).toBe('');
        const notice = annotationPanel.getByRole('status');
        await expect(notice).toContainText('notes, gene_desc, pmid');
        await expect(notice).toContainText('not saved in Sessions or TSV re-export');
        await expect(notice).toHaveAttribute('aria-live', 'polite');
        await notice.evaluate(element => element.scrollIntoView({ block: 'center' }));
        expect(await notice.evaluate(element => element.scrollWidth <= element.clientWidth)).toBe(true);
        await page.screenshot({ path: testInfo.outputPath('auxiliary-notice.png') });
        expect((await inspect(page, mode)).runs).toBe(source.runs);
        const colorInput = page.getByLabel('Specific Table (-t)', { exact: true });
        await reveal(colorInput);
        const beforeImport = await inspect(page, mode);
        const originalColorFile = await colorFile(page);
        await colorInput.setInputFiles({ name: 'integration-colors.tsv', mimeType: 'text/plain', buffer: Buffer.from(input.colors) });
        await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
        const imported = await inspect(page, mode);
        expect(imported.undoCount).toBe(beforeImport.undoCount + 1);
        await history(page, 'Undo');
        const importUndone = await inspect(page, mode);
        expect(await drawingSignature(page, importUndone.svg)).toBe(await drawingSignature(page, beforeImport.svg));
        expect(importUndone.rules).toEqual(beforeImport.rules);
        expect(importUndone.warnings).toEqual(beforeImport.warnings);
        expect(importUndone.request).toEqual(beforeImport.request);
        expect(importUndone.geometry).toEqual(beforeImport.geometry);
        expect(await colorFile(page)).toEqual(originalColorFile);
        await history(page, 'Redo');
        const importRedone = await inspect(page, mode);
        expect(await drawingSignature(page, importRedone.svg)).toBe(await drawingSignature(page, imported.svg));
        expect(importRedone.rules).toEqual(imported.rules);
        expect(await colorFile(page)).toEqual({ name: 'integration-colors.tsv', text: input.colors });
        if (mode === 'linear') await page.evaluate(setId => {
          const slots = window.__GBDRAW_APP__.adv.linear_track_slots;
          slots.push({ id: 'integration_annotations', renderer: 'annotations', enabled: true, side: 'above', height: '40px', spacing: '5px', z: 0, params: { set_id: setId, overflow: 'error', show_labels: true } });
        }, input.setId);
        const panel = await openStack(page, mode);
        await panel.locator('..').getByLabel('Use custom stack', { exact: true }).check();
        const row = page.locator(`[data-capture="${mode}-track-slot-gc_content"]`);
        const titles = mode === 'linear' ? ['Height', 'Spacing'] : ['Inner gap', 'Outer gap'];
        for (const title of titles) { const field = row.getByTitle(title, { exact: true }); await field.fill('10 PX'); await field.press('Tab'); }
        await generate(page);
        console.log(`integration ${mode}/${width}: combined generated`);
        const combined = await inspect(page, mode);
        expect(combined.features).toBe(source.features);
        expect(combined.requestSent).toEqual(combined.request);
        expect(combined.warnings).toHaveLength(1);
        expect(combined.warnings[0]).toMatchObject({ code: 'feature_selector_unmatched', annotationId: 'integration-miss', missingCount: 1 });
        expect(combined.svg).not.toContain('SKIPPED-INTEGRATION');
        expect(combined.svg).toContain('data-gbdraw-annotation-id');
        expect(JSON.stringify(combined.request)).not.toContain('PRIVATE-AUX-CELL');
        const slot = combined.request.diagramOptions.tracks[`${mode}TrackSlots`].find(slot => slot.id === 'gc_content');
        if (mode === 'linear') { expect(slot.height).toEqual({ value: 10, unit: 'px' }); expect(slot.spacing).toEqual({ value: 10, unit: 'px' }); }
        else { expect(slot.innerGapPx).toBe(10); expect(slot.outerGapPx).toBe(10); }
        for (const legend of combined.legends) {
          for (const color of ['#112233', '#445566']) expect(legend).toContainEqual({ caption: `Integration group [${color}]`, color });
        }
        const status = page.getByTestId('annotation-resolution-notice');
        await expect(status).toContainText(`${input.setId}/integration-miss`);
        await status.evaluate(element => element.scrollIntoView({ block: 'center' }));
        expect(await status.evaluate(element => element.scrollWidth <= element.clientWidth)).toBe(true);
        await page.screenshot({ path: testInfo.outputPath('combined-warning.png') });
        // Canonical live edit and its artifact warnings travel in the same History transaction.
        await page.evaluate(async id => {
          const a = window.__GBDRAW_APP__;
          await a.setSpecificRuleField(a.manualSpecificRules.findIndex(rule => rule.val === id), 'cap', 'Integration primary');
        }, input.ids[0]);
        const live = await inspect(page, mode);
        expect(live.undoCount).toBe(combined.undoCount + 1);
        expect(live.legends[0]).toContainEqual({ caption: 'Integration primary', color: '#112233' });
        expect(live.warnings).toEqual(combined.warnings);
        await history(page, 'Undo');
        const undone = await inspect(page, mode);
        await fs.writeFile(testInfo.outputPath('before-undo.svg'), combined.svg);
        await fs.writeFile(testInfo.outputPath('after-undo.svg'), undone.svg);
        await fs.writeFile(testInfo.outputPath('undo-state.json'), JSON.stringify({ before: { ...combined, svg: undefined }, after: { ...undone, svg: undefined } }, null, 2));
        expect(await drawingSignature(page, undone.svg)).toBe(await drawingSignature(page, combined.svg));
        expect(undone.rules).toEqual(combined.rules);
        expect(undone.legends).toEqual(combined.legends);
        expect(undone.warnings).toEqual(combined.warnings);
        expect(undone.request).toEqual(combined.request);
        expect(undone.geometry).toEqual(combined.geometry);
        await history(page, 'Redo');
        const redone = await inspect(page, mode);
        await fs.writeFile(testInfo.outputPath('live.svg'), live.svg);
        await fs.writeFile(testInfo.outputPath('redone.svg'), redone.svg);
        await fs.writeFile(testInfo.outputPath('redo-state.json'), JSON.stringify({ live: { ...live, svg: undefined }, redone: { ...redone, svg: undefined } }, null, 2));
        expect(await drawingSignature(page, redone.svg)).toBe(await drawingSignature(page, live.svg));
        expect(redone.rules).toEqual(live.rules);
        expect(redone.legends).toEqual(live.legends);
        expect(redone.warnings).toEqual(live.warnings);
        expect(redone.request).toEqual(live.request);
        expect(redone.geometry).toEqual(live.geometry);
        // Invalid pixels do not coerce a combined Result or its warning to a new artifact.
        const field = row.getByTitle(titles[0], { exact: true });
        await field.fill('10em'); await field.press('Tab');
        await expect(row.getByRole('alert')).toContainText('finite');
        const beforeFailure = await inspect(page, mode);
        await page.getByRole('button', { name: 'Generate Diagram', exact: true }).press('Enter');
        await expect.poll(async () => { const current = await inspect(page, mode); return { completed: current.completed, processing: current.processing }; }, { timeout: 180000 })
          .toEqual({ completed: beforeFailure.completed + 1, processing: false });
        const failed = await inspect(page, mode);
        expect(artifact(failed)).toEqual(artifact(redone)); expect(failed.runs).toBe(beforeFailure.runs);
        expect(failed.slots.find(slot => slot.id === 'gc_content')[mode === 'linear' ? 'height' : 'inner_gap_px']).toBe('10em');
        await row.getByRole('alert').evaluate(element => element.scrollIntoView({ block: 'center' }));
        await page.screenshot({ path: testInfo.outputPath('combined-pixel-error.png') });
        await field.fill('1e1px');
        await generate(page);
        const final = await inspect(page, mode);
        await row.getByTitle('Enabled', { exact: true }).uncheck();
        const sessionPath = testInfo.outputPath('combined-disabled.gbdraw-session.json.gz');
        const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', sessionPath)));
        expect(saved.runMetadata.annotationWarnings).toEqual(final.warnings);
        expect(saved.renderRequest).toEqual(final.request);
        expect(saved.config.adv[`${mode}_track_slots`].find(slot => slot.id === 'gc_content').enabled).toBe(false);
        expect(JSON.stringify(saved.config.annotationSets)).not.toContain('PRIVATE-AUX-CELL');
        restored = await load(browser, sessionPath, { width, height: 1000 });
        const loaded = await inspect(restored, mode);
        expect(artifact(loaded)).toEqual(artifact(final)); expect(loaded.rules).toEqual(final.rules);
        expect((await getDiagramWorkerActivity(restored)).runs).toBe(0);
        await openStack(restored, mode);
        await restored.locator(`[data-capture="${mode}-track-slot-gc_content"]`).getByTitle('Disabled', { exact: true }).check();
        await generate(restored);
        const regenerated = await inspect(restored, mode);
        expect(artifact(regenerated)).toEqual(artifact(final));
        expect(await restored.evaluate(async () => { const { readFileText } = await import('./js/services/file-content-cache.js'); return readFileText(window.__GBDRAW_APP__.files.t_color); })).toBe(input.colors);
        const activePath = testInfo.outputPath('combined.gbdraw-session.json.gz');
        const active = JSON.parse(gunzipSync(await download(restored, 'Save Session', activePath)));
        expect(active.resources).toEqual(saved.resources);
        const svgPath = testInfo.outputPath('combined.svg');
        const pendingSvg = restored.waitForEvent('download');
        await restored.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
        await (await pendingSvg).saveAs(svgPath);
        expect(await fs.readFile(svgPath, 'utf8')).not.toContain('SKIPPED-INTEGRATION');
        const pendingRules = restored.waitForEvent('download');
        await restored.evaluate(() => window.__GBDRAW_APP__.downloadSpecificRulesTsv());
        await (await pendingRules).saveAs(testInfo.outputPath('combined-rules.tsv'));
        expect(await fs.readFile(testInfo.outputPath('combined-rules.tsv'), 'utf8')).toContain('\tIntegration primary\n');
        const pendingAnnotations = page.waitForEvent('download');
        await annotationPanel.getByRole('button', { name: 'Download TSV', exact: true }).click();
        const annotationPath = testInfo.outputPath('combined-annotations.tsv');
        await (await pendingAnnotations).saveAs(annotationPath);
        const annotationText = await fs.readFile(annotationPath, 'utf8');
        expect(annotationText).toContain('PRIVATE-INTEGRATION-MISSING');
        expect(annotationText).not.toContain('PRIVATE-AUX-CELL');
        expect(annotationText.split('\n')[0].split('\t')).not.toContain('notes');
        const native = JSON.parse(execFileSync('python', ['-c', `
import json,sys,xml.etree.ElementTree as E
from pathlib import Path
from gbdraw.session import materialize_session,render_session
from gbdraw.render.track_slot_metadata import collect_track_slot_geometry_records,build_track_slot_geometry_run_metadata
with materialize_session(sys.argv[1],output_directory=sys.argv[2]) as session:
 result=render_session(session)
 svg=result.drawing.tostring()
 Path(sys.argv[2],'native.svg').write_text(svg)
 root=E.fromstring(svg)
 captions={n.get('data-legend-key'):next((p.get('fill') for p in n.iter() if p.tag.endswith('path') and p.get('fill') not in (None,'none')),None) for n in root.iter() if n.get('data-legend-key')}
 geometry=build_track_slot_geometry_run_metadata(mode=result.mode,records=collect_track_slot_geometry_records(result.drawing,result_index=0,result_name=result.output_paths[0].name))['trackSlotGeometry']
 print(json.dumps({'geometry':geometry,'captions':captions,'rules':result.request.options.colors.color_table['caption'].tolist(),'warnings':[w.code for w in result.annotation_warnings],'skippedMark':'SKIPPED-INTEGRATION' in svg}))
`, activePath, testInfo.outputPath('native')], { encoding: 'utf8', maxBuffer: 4 * 1024 * 1024 }));
        expect(native.geometry).toEqual(final.geometry);
        expect(native.rules).toEqual(final.rules.map(rule => rule.cap));
        expect(native.warnings).toEqual(['feature_selector_unmatched']); expect(native.skippedMark).toBe(false);
        for (const entry of regenerated.legends[0].filter(entry => entry.caption.startsWith('Integration'))) expect(native.captions[entry.caption]).toBe(entry.color);
        expect(consoleText.join('\n')).not.toMatch(/PRIVATE-INTEGRATION-MISSING|PRIVATE-AUX-CELL|PRIVATE-DESC/);
        expect(page.externalRequests).toEqual([]); expect(restored.externalRequests).toEqual([]);
        await fs.writeFile(testInfo.outputPath('integration-evidence.json'), JSON.stringify({ seed, width, mode, sourceFeatures: source.features,
          request: final.request, warnings: final.warnings, geometry: final.geometry, native, version: active.version,
          originalColorBytesPreserved: true, lazyLoadRenderRuns: 0, externalRequests: page.externalRequests }, null, 2));
      } finally { await page.context().close(); if (restored) await restored.context().close(); }
    });
  }
}
