const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { execFileSync } = require('node:child_process');
const { load, generate, download, loadEditorLegendRows } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

// A rule edit that changes the Legend asks for the automatic rerender (OV-42,
// OV-43, #857), so the Result is read once Python has drawn it.
const inspect = async page => {
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing && !window.__GBDRAW_APP__.labelReflowProcessing);
  await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing && !window.__GBDRAW_APP__.labelReflowProcessing);
  return inspectNow(page);
};
const inspectNow = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const { getVisibleFeatureLegendGroup, getAllFeatureLegendGroups, legendRowShown } = await import('./js/services/legend-svg.js');
  const root=s.svgContainer.value.querySelector('svg');
  const result=new DOMParser().parseFromString(s.results.value[s.selectedResultIndex.value].content,'image/svg+xml').documentElement;
  // The rows a Result shows: a rule commit that merges a row into a shown row
  // retires it live as a hidden row, which the layout skips and every export
  // strips, until Python draws the Legend again (D-15-6 (2), (3); OV-313).
  const entries=svg=>[...(getVisibleFeatureLegendGroup(svg)?.querySelectorAll('g[data-legend-key]')||[])].filter(legendRowShown).map(e=>({caption:e.getAttribute('data-legend-key'),color:e.querySelector('path[fill]')?.getAttribute('fill')}));
  // The name each row draws, by its key (`caption` above is the key, which is never drawn).
  const names=svg=>Object.fromEntries([...(getVisibleFeatureLegendGroup(svg)?.querySelectorAll('g[data-legend-key]')||[])].filter(legendRowShown).map(e=>[e.getAttribute('data-legend-key'),e.querySelector('text')?.textContent]));
  return {rules:JSON.parse(JSON.stringify(s.activeDrawing().manualSpecificRules.map(rule=>({...rule,fromFile:Boolean(rule.fromFile)})))), mounted:entries(root), result:entries(result),
    dual:getAllFeatureLegendGroups(root).map(group=>[...group.querySelectorAll('g[data-legend-key]')].map(e=>e.getAttribute('data-legend-key'))),
    dualStyles:getAllFeatureLegendGroups(root).map(group=>[...group.querySelectorAll('g[data-legend-key]')].map(e=>({caption:e.getAttribute('data-legend-key'),color:e.querySelector('path[fill]')?.getAttribute('fill')}))),
    names:names(result), mountedNames:names(root),
    svg:s.results.value[s.selectedResultIndex.value].content, fileName:s.files.t_color?.name};
});
const reveal = async locator => {
  for(const details of await locator.locator('xpath=ancestor::details').all()) {
    if(await details.getAttribute('open')===null) await details.locator(':scope > summary').press('Enter');
  }
};
const history = async (page, name) => {
  await page.getByRole('button',{name,exact:true}).click();
  await expect.poll(()=>page.evaluate(()=>!window.__GBDRAW_HISTORY__.restoring.value && !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
};

for(const width of [1440,390]) {
  test(`a caption written in several colors stays as written and draws numbered names through live/native/edit/session/download (${width}px)`,async({browser},testInfo)=>{
    test.setTimeout(1_800_000);
    const page=await load(browser,'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json',{width,height:1000});
    let fresh;
    try {
      page.setDefaultTimeout(180000);
      await generate(page);
      const source=await evaluateWithRetainedPromise(page, async()=>{
        const a=window.__GBDRAW_APP__;
        const {getFeatureColorRuleHash}=await import('./js/services/feature-utils.js');
        const selected=a.extractedFeatures.filter(f=>f.type==='CDS').slice(0,3);
        const ids=selected.map(getFeatureColorRuleHash);
        const original=a.manualSpecificRules.map(r=>[r.feat,r.qual,r.val,r.color,r.cap].join('\t')).join('\n');
        Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:ids[2],color:'#778899',cap:'Color group'});
        await a.addSpecificRule();
        return {ids,text:original+'\n'+ids.slice(0,2).map((id,i)=>['CDS','hash',id,['#112233','#445566'][i],'Color group'].join('\t')).join('\n')+'\n'};
      });
      const before=await inspect(page);
      const input=page.getByLabel('Specific rules file',{exact:true});
      await reveal(input);
      const picker=page.getByRole('button',{name:'Choose Specific rules file',exact:true});
      await picker.scrollIntoViewIfNeeded();
      await expect(picker).toBeVisible();
      const chooser=page.waitForEvent('filechooser');
      await picker.press('Enter');
      await (await chooser).setFiles({name:'multicolor.tsv',mimeType:'text/plain',buffer:Buffer.from(source.text)});
      await page.evaluate(async()=>{const a=window.__GBDRAW_APP__;await a.waitForAuxiliaryFileImport(a.files.t_color);});
      const live=await inspect(page);
      // PD-OI-042 scenario revision 2: the rules keep their caption; each color has
      // its own row (the keys below, never drawn), and the Legend numbers the names.
      const captions=['Color group [#112233]','Color group [#445566]','Color group [#778899]'];
      expect(live.rules.filter(r=>r.cap.startsWith('Color group')).map(r=>r.cap)).toEqual(['Color group','Color group','Color group']);
      expect(captions.map(key=>live.names[key]).sort()).toEqual(['Color group (1)','Color group (2)','Color group (3)']);
      expect(live.mountedNames).toEqual(live.names);
      expect(Object.values(live.names).filter(name=>name.includes(' [#'))).toEqual([]);
      expect(live.rules.find(r=>r.val===source.ids[2]).fromFile).toBe(false);
      expect(live.rules.find(r=>r.val===source.ids[0]).fromFile).toBe(true);
      expect(live.mounted).toEqual(live.result);
      for(const caption of captions) expect(live.result.some(e=>e.caption===caption)).toBe(true);
      expect(live.dual).toHaveLength(1);
      const fills=await page.evaluate(ids=>{
        const a=window.__GBDRAW_APP__,svg=a.svgContainer.querySelector('svg');
        return ids.map(id=>[...new Set([...svg.querySelectorAll(`[data-gbdraw-feature-id="${id}"]`)]
          .filter(element=>element.hasAttribute('fill')).map(element=>element.getAttribute('fill')))]);
      },source.ids);
      expect(fills).toEqual([['#112233'],['#445566'],['#778899']]);
      const extents=await page.evaluate(captions=>{
        const svg=window.__GBDRAW_APP__.svgContainer.querySelector('svg'),bounds=svg.getBoundingClientRect();
        return captions.map(caption=>{
          const entry=svg.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
          const box=entry.getBoundingClientRect(),swatch=entry.querySelector('path[fill]').getBBox();
          return {width:box.width,height:box.height,swatchWidth:swatch.width,swatchHeight:swatch.height,
            inside:box.x>=bounds.x-1 && box.y>=bounds.y-1 && box.right<=bounds.right+1 && box.bottom<=bounds.bottom+1};
        });
      },captions);
      for(const extent of extents){expect(extent.inside).toBe(true);expect(extent.width).toBeGreaterThan(0);
        expect(extent.height).toBeGreaterThan(0);expect(extent.swatchWidth).toBeGreaterThan(0);expect(extent.swatchHeight).toBeGreaterThan(0);}
      // No caption is rewritten, so nothing is announced.
      await expect(page.getByText('specific-color caption')).toHaveCount(0);
      await history(page,'Undo');
      expect((await inspect(page)).rules).toEqual(before.rules);
      await history(page,'Redo');
      expect((await inspect(page)).rules).toEqual(live.rules);
      await generate(page);
      const generated=await inspect(page);
      expect(captions.map(key=>generated.names[key])).toEqual(captions.map(key=>live.names[key]));
      expect(generated.result.filter(e=>e.caption.startsWith('Color group')).sort((a,b)=>a.caption.localeCompare(b.caption))).toEqual(live.result.filter(e=>e.caption.startsWith('Color group')).sort((a,b)=>a.caption.localeCompare(b.caption)));
      // Rename/recolor acts on the canonical source rows and leaves siblings intact.
      await page.evaluate(async id=>{
        const a=window.__GBDRAW_APP__;
        const index=a.manualSpecificRules.findIndex(r=>r.val===id);
        await a.setSpecificRuleField(index,'cap','Primary complex');
        await a.updateLegendEntryColor(a.legendEntries.findIndex(e=>e.caption==='Primary complex'),'#224466');
        a.sortLegendEntries();
      },source.ids[0]);
      const edited=await inspect(page);
      expect(edited.rules.find(r=>r.val===source.ids[0]).cap).toBe('Primary complex');
      expect(edited.rules.find(r=>r.val===source.ids[0]).color).toBe('#224466');
      expect(edited.rules.find(r=>r.val===source.ids[1]).cap).toBe('Color group');
      await generate(page);
      const final=await inspect(page);
      expect(final.result.find(e=>e.caption==='Primary complex').color).toBe('#224466');
      // Two colors are left under the name: their rows are numbered from 1.
      expect(captions.slice(1).map(key=>final.names[key])).toEqual(['Color group (1)','Color group (2)']);
      const pendingTsv=page.waitForEvent('download');
      await page.evaluate(()=>window.__GBDRAW_APP__.downloadSpecificRulesTsv());
      const tsvPath=testInfo.outputPath('canonical-rules.tsv');
      await (await pendingTsv).saveAs(tsvPath);
      const tsv=await fs.readFile(tsvPath,'utf8');
      expect(tsv).toContain('\tPrimary complex\n');
      expect(tsv).toContain('\t#445566\tColor group\n');
      expect(tsv).not.toContain(' [#');
      expect(await page.evaluate(async()=>{const {readFileText}=await import('./js/services/file-content-cache.js');const {state}=await import('./js/state.js');return readFileText(state.inputFiles.fileOf(window.__GBDRAW_APP__.files.t_color));})).toBe(source.text);
      expect(final.fileName).toBe('multicolor.tsv');
      const sessionPath=testInfo.outputPath('canonical.gbdraw-session.json.gz');
      await download(page,'Save Session',sessionPath);
      const svgPath=testInfo.outputPath('canonical.svg');
      const svgPending=page.waitForEvent('download');
      await page.evaluate(()=>window.__GBDRAW_APP__.downloadSVG());
      await (await svgPending).saveAs(svgPath);
      fresh=await load(browser,sessionPath,{width,height:1000});
      const loaded=await inspect(fresh);
      expect(loaded.rules).toEqual(final.rules);
      expect(loaded.result).toEqual(final.result);
      expect((await getDiagramWorkerActivity(fresh)).runs).toBe(0);
      await generate(fresh);
      expect((await inspect(fresh)).result.filter(e=>e.caption==='Primary complex'||e.caption.startsWith('Color group'))).toEqual(final.result.filter(e=>e.caption==='Primary complex'||e.caption.startsWith('Color group')));
      const native=JSON.parse(execFileSync('python',['-c',`
import json,sys,xml.etree.ElementTree as E
from gbdraw.session import materialize_session,render_session
with materialize_session(sys.argv[1],output_directory=sys.argv[2]) as session:
 result=render_session(session)
 # The Gallery Session keeps its command line's grouping: one Circular diagram per record.
 [drawing]=getattr(result,'drawings',None) or (result.drawing,)
 root=E.fromstring(drawing.tostring())
 captions={n.get('data-legend-key'):next((p.get('fill') for p in n.iter() if p.tag.endswith('path') and p.get('fill') not in (None,'none')),None) for n in root.iter() if n.get('data-legend-key')}
 print(json.dumps({'captions':captions,'rules':result.request.options.colors.color_table['caption'].tolist()}))
`,sessionPath,testInfo.outputPath('native')],{encoding:'utf8',maxBuffer:4*1024*1024}));
      expect(native.rules).toEqual(final.rules.map(r=>r.cap));
      for(const entry of final.result.filter(e=>e.caption==='Primary complex'||e.caption.startsWith('Color group'))) expect(native.captions[entry.caption]).toBe(entry.color);
      await fs.writeFile(testInfo.outputPath('caption-evidence.json'),JSON.stringify({source,live: {...live,svg:undefined},final:{...final,svg:undefined},native},null,2));
      await input.setInputFiles({name:'multicolor-reimport.tsv',mimeType:'text/plain',buffer:Buffer.from(source.text)});
      await page.evaluate(async()=>{const a=window.__GBDRAW_APP__;await a.waitForAuxiliaryFileImport(a.files.t_color);});
      const reimported=await inspect(page);
      expect(reimported.rules.find(r=>r.val===source.ids[0] && !r.fromFile).cap).toBe('Primary complex');
      expect(reimported.result.some(e=>e.caption==='Color group [#112233]')).toBe(false);
      await page.evaluate(async()=>{const a=window.__GBDRAW_APP__;await a.removeSpecificRule(a.manualSpecificRules.findIndex(r=>r.cap==='Color group' && r.color==='#445566'));});
      expect((await inspect(page)).result.some(e=>e.caption==='Color group [#445566]')).toBe(false);
      await history(page,'Undo');
      expect((await inspect(page)).rules).toEqual(reimported.rules);
      await history(page,'Redo');
      expect((await inspect(page)).result.some(e=>e.caption==='Color group [#445566]')).toBe(false);
      expect((await inspect(page)).result.find(e=>e.caption==='Primary complex').color).toBe('#224466');
      expect(page.externalRequests).toEqual([]);
    }finally{await page.context().close();if(fresh)await fresh.context().close();}
  });
}


test('a rule row whose key is the key of another listed row is refused, and a Session draft of one caption in two colors keeps it',async({browser},testInfo)=>{
  test.setTimeout(1_800_000);
  const page=await load(browser,'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
  page.setDefaultTimeout(180000);
  let fresh;
  try {
    await generate(page);
    // An editor row a Session holds (R15-2 retired Add legend item), named like
    // the key of the row the second rule below would draw.
    await loadEditorLegendRows(page,[['Conflict [#112233]','#abcdef']]);
    await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.legendEntries.some(e=>e.caption==='Conflict [#112233]'))).toBe(true);
    const ids=await evaluateWithRetainedPromise(page, async()=>{
      const a=window.__GBDRAW_APP__,{getFeatureColorRuleHash}=await import('./js/services/feature-utils.js');
      const ids=a.extractedFeatures.filter(f=>f.type==='CDS').slice(0,2).map(getFeatureColorRuleHash);
      Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:ids[0],color:'#445566',cap:'Conflict'});
      await a.addSpecificRule();return ids;
    });
    const before=await inspect(page);
    const errorLog=()=>page.evaluate(async()=>{const e=(await import('./js/state.js')).state.errorLog.value;return e&&{code:e.code,operation:e.operation,summary:e.summary};});
    const warningBefore=await page.evaluate(async()=>JSON.stringify((await import('./js/state.js')).state.annotationWarnings.value));
    expect(await errorLog()).toBeNull();
    // D-16: a second color of the caption gives each color its own row key
    // ("Conflict [#112233]"), which the editor row above already shows. No two
    // rows share one shown key: the commit is refused and nothing is written.
    await evaluateWithRetainedPromise(page, async id=>{
      const a=window.__GBDRAW_APP__;
      Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:id,color:'#112233',cap:'Conflict'});
      await a.addSpecificRule();
    },ids[1]);
    const rejected=await inspect(page);
    await fs.writeFile(testInfo.outputPath('collision-evidence.json'),JSON.stringify({before:{...before,svg:undefined},rejected:{...rejected,svg:undefined},error:await errorLog()},null,2));
    expect(rejected.rules).toEqual(before.rules);expect(rejected.svg).toBe(before.svg);
    expect(rejected.result).toEqual(before.result);expect(rejected.names).toEqual(before.names);
    expect(rejected.names['Conflict']).toBe('Conflict');expect(rejected.names['Conflict [#112233]']).toBe('Conflict [#112233]');
    // OV-644: the refusal says what to change, without the row key.
    expect(await errorLog()).toEqual({code:'LEGEND_ROW_CONFLICT',operation:'evaluateRules',
      summary:'A color rule would draw a Legend item that already exists in another color. Change the rule\'s Legend name or color, or rename that Legend item.'});
    expect(await page.evaluate(async()=>JSON.stringify((await import('./js/state.js')).state.annotationWarnings.value))).toBe(warningBefore);
    // A draft fixture: one listed row of a caption that two rules write in two
    // colors, with an unchanged preview.
    // The row is renamed before its name is a color rule's: a rename onto the
    // name of a color rule asks Merge / Keep both / Cancel (PD-OI-061 revision 4).
    await page.evaluate(async ids=>{
      const a=window.__GBDRAW_APP__;
      const {manualSpecificRules,featureColorOverrides}=a;
      await a.renameLegendEntry(a.legendEntries.findIndex(entry=>entry.caption==='Conflict'),'Historical');
      manualSpecificRules.splice(0,manualSpecificRules.length,
        ...ids.map((id,i)=>({feat:'CDS',qual:'hash',val:id,color:['#112233','#445566'][i],cap:'Historical'})));
      for(const override of Object.values(featureColorOverrides)) if(override.caption==='Conflict') override.caption='Historical';
      const {state:s}=await import('./js/state.js'),{serializeCleanSvg}=await import('./js/services/svg-serialization.js');
      const svg=s.svgContainer.value.querySelector('svg');
      for(const entry of svg.querySelectorAll('g[data-legend-key="Historical"]')) entry.removeAttribute('data-legend-owner');
      const results=[...s.results.value];results[s.selectedResultIndex.value]={...results[s.selectedResultIndex.value],content:serializeCleanSvg(svg)};
      s.results.value=results;
    },ids);
    const draft=await inspect(page),saved=testInfo.outputPath('historical.gbdraw-session.json.gz');
    await download(page,'Save Session',saved);
    fresh=await load(browser,saved);fresh.setDefaultTimeout(180000);
    const loaded=await inspect(fresh);
    expect(loaded.rules.map(r=>r.cap)).toEqual(['Historical','Historical']);
    expect(loaded.result).toEqual(draft.result);expect((await getDiagramWorkerActivity(fresh)).runs).toBe(0);
    await fresh.evaluate(async id=>{await window.__GBDRAW_APP__.setSpecificRuleField(0,'val',id);},ids[0]);
    const edited=await inspect(fresh);
    await fs.writeFile(testInfo.outputPath('draft-evidence.json'),JSON.stringify({draft:{...draft,svg:undefined},loaded:{...loaded,svg:undefined},edited:{...edited,svg:undefined}},null,2));
    // PD-OI-042 scenario revision 2: the captions stay as written; each color
    // has its own row (the keys below, never drawn), numbered in the Legend.
    const keys=[['Historical [#112233]','#112233'],['Historical [#445566]','#445566']];
    expect(edited.rules.map(r=>r.cap)).toEqual(['Historical','Historical']);
    for(const [key,color] of keys) expect(edited.result.find(e=>e.caption===key)?.color).toBe(color);
    expect(edited.result.some(e=>e.caption==='Historical')).toBe(false);
    expect(keys.map(([key])=>edited.names[key]).sort()).toEqual(['Historical (1)','Historical (2)']);
    expect(edited.mountedNames).toEqual(edited.names);
    expect(edited.result.find(e=>e.caption==='Conflict [#112233]').color).toBe('#abcdef');
    expect(edited.names['Conflict [#112233]']).toBe('Conflict [#112233]');
    await generate(fresh);
    const generated=await inspect(fresh);
    expect(generated.rules.map(r=>r.cap)).toEqual(['Historical','Historical']);
    expect(generated.result.filter(e=>e.caption.startsWith('Historical'))).toEqual(edited.result.filter(e=>e.caption.startsWith('Historical')));
    expect(keys.map(([key])=>generated.names[key])).toEqual(keys.map(([key])=>edited.names[key]));
  } finally {await page.context().close();if(fresh)await fresh.context().close();}
});


test('Linear comparison numbers a caption written in two colors in both legend orientations',async({browser},testInfo)=>{
  test.setTimeout(1_800_000);
  const page=await load(browser,'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
  page.setDefaultTimeout(180000);
  try {
    await generate(page);
    await evaluateWithRetainedPromise(page, async()=>{
      const a=window.__GBDRAW_APP__,{getFeatureColorRuleHash}=await import('./js/services/feature-utils.js');
      const features=a.extractedFeatures.filter(f=>f.type==='CDS').slice(0,2);
      for(let i=0;i<features.length;i++) {
        Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:getFeatureColorRuleHash(features[i]),color:['#112233','#445566'][i],cap:'Comparison group'});
        await a.addSpecificRule();
      }
    });
    const live=await inspect(page);
    expect(live.dual).toHaveLength(2);expect([...live.dual[0]].sort()).toEqual([...live.dual[1]].sort());
    for(const caption of ['Comparison group [#112233]','Comparison group [#445566]'])
      for(const orientation of live.dual) expect(orientation).toContain(caption);
    for(const orientation of live.dualStyles) expect(orientation.filter(e=>e.caption.startsWith('Comparison group'))).toEqual([
      {caption:'Comparison group [#112233]',color:'#112233'}, {caption:'Comparison group [#445566]',color:'#445566'}]);
    expect(live.mounted).toEqual(live.result);
    // The rules keep the caption; both orientations draw the numbered names.
    expect(live.rules.filter(r=>r.cap.startsWith('Comparison group')).map(r=>r.cap)).toEqual(['Comparison group','Comparison group']);
    const drawn=await page.evaluate(async()=>{
      const {getAllFeatureLegendGroups}=await import('./js/services/legend-svg.js');
      return getAllFeatureLegendGroups(window.__GBDRAW_APP__.svgContainer.querySelector('svg')).map(group=>
        [...group.querySelectorAll('g[data-legend-key^="Comparison group"]')].map(e=>e.querySelector('text')?.textContent));
    });
    expect(drawn).toEqual([['Comparison group (1)','Comparison group (2)'],['Comparison group (1)','Comparison group (2)']]);
    await generate(page);
    const fresh=await inspect(page);
    expect(fresh.result.filter(e=>e.caption.startsWith('Comparison group'))).toEqual(live.result.filter(e=>e.caption.startsWith('Comparison group')));
    expect(['Comparison group [#112233]','Comparison group [#445566]'].map(key=>fresh.names[key])).toEqual(['Comparison group (1)','Comparison group (2)']);
    const pending=page.waitForEvent('download');
    await page.evaluate(()=>window.__GBDRAW_APP__.downloadSVG());
    await (await pending).saveAs(testInfo.outputPath('comparison.svg'));
    await page.locator('.origin-top svg').screenshot({path:testInfo.outputPath('comparison.svg.png')});
    await fs.writeFile(testInfo.outputPath('dual-legend.json'),JSON.stringify({...fresh,svg:undefined},null,2));
  }finally{await page.context().close();}
});


for (const width of [1440, 390]) {
  // D-15: "Apply to all" on the loaded tRNA palette row sets the tRNA default
  // color, so the whole row, hidden features included, takes it at once.
  test(`loaded palette row takes its new default color atomically (${width}px)`, async ({ browser }, info) => {
    test.setTimeout(180000);
    const page = await load(browser, 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', { width, height: 1000 });
    let fresh;
    try {
      const evidence = () => page.evaluate(async () => {
        const a = window.__GBDRAW_APP__;
        const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
        const root = a.svgContainer.querySelector('svg');
        const result = new DOMParser().parseFromString(a.results[a.selectedResultIndex].content, 'image/svg+xml');
        const ids = new Set(a.extractedFeatures.filter(f => f.type === 'tRNA').map(f => f.svg_id));
        const colors = svg => [...new Set([...svg.querySelectorAll('path[data-gbdraw-feature-part="block"]')].filter(e => ids.has(e.getAttribute('data-gbdraw-feature-id'))).map(e => e.getAttribute('fill')))];
        return { mounted: colors(root), result: colors(result),
          undo: window.__GBDRAW_HISTORY__.getUndoCount(),
          committedScaleInterval: getCommittedCanonicalRenderRequest()?.diagramOptions?.configOverrides?.['objects.scale.interval'] ?? null };
      });
      const target = await page.evaluate(() => {
        const a = window.__GBDRAW_APP__, feature = a.extractedFeatures.find(f => f.type === 'tRNA');
        a.openFeatureEditorFromList(feature, { clientX: 220, clientY: 220 });
        a.adv.scale_interval = 12345;
        return feature.svg_id;
      });
      expect(target).toBeTruthy();
      const before = await evidence();
      expect(before.mounted).toHaveLength(1);
      expect(before.result).toEqual(before.mounted);
      const fill = page.getByRole('dialog', { name: /Feature details:/ }).getByLabel('Feature fill color').first();
      await fill.fill('#c026d3');
      const scope = page.getByRole('heading', { name: 'Color Change Scope' }).locator('..');
      await scope.getByRole('button').filter({ hasText: 'Legend item "tRNA"' }).last().click();
      // The choice prepares the rules it adds; on a Session just loaded that
      // starts the Worker runtime (OV-225), so wait for the choice to finish.
      await page.waitForFunction(() => !window.__GBDRAW_APP__.featureStyleScopeDialog.show);
      await expect.poll(async () => (await evidence()).mounted, { timeout: 10000 }).toEqual(['#c026d3']);
      const after = await evidence();
      expect(after.result).toEqual(['#c026d3']);
      expect(after.undo).toBe(before.undo + 1);
      // The scale draft stays out of the committed request; no derived status exists.
      expect(after.committedScaleInterval).toBe(before.committedScaleInterval);
      expect(after.committedScaleInterval).not.toBe(12345);
      await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);
      const live = await inspect(page);
      expect(live.result.find(entry => entry.caption === 'tRNA')?.color).toBe('#c026d3');
      expect(live.rules.filter(rule => rule.feat === 'tRNA')).toEqual([]);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.currentColors.tRNA)).toBe('#c026d3');
      await history(page, 'Undo');
      expect((await evidence()).mounted).toEqual(before.mounted);
      expect((await evidence()).result).toEqual(before.result);
      expect((await inspect(page)).result.find(entry => entry.caption === 'tRNA')?.color).toBe(before.mounted[0]);
      await history(page, 'Redo');
      expect((await evidence()).mounted).toEqual(after.mounted);
      expect((await inspect(page)).rules).toEqual(live.rules);
      expect((await inspect(page)).result.find(entry => entry.caption === 'tRNA')?.color).toBe('#c026d3');
      const saved = info.outputPath('default-caption.json.gz');
      await download(page, 'Save Session', saved);
      fresh = await load(browser, saved, { width, height: 1000 });
      expect((await inspect(fresh)).rules).toEqual(live.rules);
      expect((await inspect(fresh)).result.find(entry => entry.caption === 'tRNA')?.color).toBe('#c026d3');
      expect((await getDiagramWorkerActivity(fresh)).runs).toBe(0);
      await generate(fresh);
      expect((await inspect(fresh)).result.find(entry => entry.caption === 'tRNA')?.color).toBe('#c026d3');
      await history(page, 'Undo');
      await evaluateWithRetainedPromise(page, async () => {
        const a = window.__GBDRAW_APP__;
        const { getFeatureColorRuleHash } = await import('./js/services/feature-utils.js');
        Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'hash',
          val: getFeatureColorRuleHash(a.extractedFeatures.find(feature => feature.type === 'CDS')),
          color: '#112233', cap: 'Independent group' });
        await a.addSpecificRule();
      });
      const independent = await inspect(page);
      expect(independent.result.find(entry => entry.caption === 'tRNA')?.color).toBe(before.mounted[0]);
      expect(independent.result.find(entry => entry.caption === 'Independent group')?.color).toBe('#112233');
      expect(page.externalRequests).toEqual([]);
      expect(fresh.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
      if (fresh) await fresh.context().close();
    }
  });
}
