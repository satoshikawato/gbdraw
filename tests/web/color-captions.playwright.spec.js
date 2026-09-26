const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { execFileSync } = require('node:child_process');
const { load, generate, download } = require('./helpers/mode-transition.cjs');
const { getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const inspect = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const { getVisibleFeatureLegendGroup, getAllFeatureLegendGroups } = await import('./js/app/legend/utils.js');
  const root=s.svgContainer.value.querySelector('svg');
  const result=new DOMParser().parseFromString(s.results.value[s.selectedResultIndex.value].content,'image/svg+xml').documentElement;
  const entries=svg=>[...(getVisibleFeatureLegendGroup(svg)?.querySelectorAll('g[data-legend-key]')||[])].map(e=>({caption:e.getAttribute('data-legend-key'),color:e.querySelector('path[fill]')?.getAttribute('fill')}));
  return {rules:JSON.parse(JSON.stringify(s.manualSpecificRules.map(rule=>({...rule,fromFile:Boolean(rule.fromFile)})))), mounted:entries(root), result:entries(result),
    dual:getAllFeatureLegendGroups(root).map(group=>[...group.querySelectorAll('g[data-legend-key]')].map(e=>e.getAttribute('data-legend-key'))),
    dualStyles:getAllFeatureLegendGroups(root).map(group=>[...group.querySelectorAll('g[data-legend-key]')].map(e=>({caption:e.getAttribute('data-legend-key'),color:e.querySelector('path[fill]')?.getAttribute('fill')}))),
    svg:s.results.value[s.selectedResultIndex.value].content, fileName:s.files.t_color?.name,
    notice:window.__GBDRAW_APP__.specificRuleNotice};
});
const reveal = async locator => {
  for(const details of await locator.locator('xpath=ancestor::details').all()) {
    if(await details.getAttribute('open')===null) await details.locator(':scope > summary').click();
  }
};
const history = async (page, name) => {
  await page.getByRole('button',{name,exact:true}).click();
  await expect.poll(()=>page.evaluate(()=>!window.__GBDRAW_HISTORY__.restoring.value && !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
};

for(const width of [1440,390]) {
  test(`canonical multicolor captions survive live/native/edit/session/download (${width}px)`,async({browser},testInfo)=>{
    test.setTimeout(1_800_000);
    const page=await load(browser,'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json',{width,height:1000});
    let fresh;
    try {
      page.setDefaultTimeout(180000);
      await generate(page);
      const source=await page.evaluate(async()=>{
        const a=window.__GBDRAW_APP__;
        const {getFeatureGenerationHash}=await import('./js/app/feature-utils.js');
        const selected=a.extractedFeatures.filter(f=>f.type==='CDS').slice(0,3);
        const ids=selected.map(getFeatureGenerationHash);
        const original=a.manualSpecificRules.map(r=>[r.feat,r.qual,r.val,r.color,r.cap].join('\t')).join('\n');
        Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:ids[2],color:'#778899',cap:'Color group'});
        await a.addSpecificRule();
        return {ids,text:original+'\n'+ids.slice(0,2).map((id,i)=>['CDS','hash',id,['#112233','#445566'][i],'Color group'].join('\t')).join('\n')+'\n'};
      });
      const before=await inspect(page);
      const input=page.getByLabel('Specific Table (-t)',{exact:true});
      await reveal(input);
      const picker=page.getByRole('button',{name:'Choose Specific Table (-t)',exact:true});
      await picker.scrollIntoViewIfNeeded();
      await expect(picker).toBeVisible();
      const chooser=page.waitForEvent('filechooser');
      await picker.press('Enter');
      await (await chooser).setFiles({name:'multicolor.tsv',mimeType:'text/plain',buffer:Buffer.from(source.text)});
      await page.evaluate(async()=>{const a=window.__GBDRAW_APP__;await a.waitForAuxiliaryFileImport(a.files.t_color);});
      const live=await inspect(page);
      const captions=['Color group [#112233]','Color group [#445566]','Color group [#778899]'];
      expect(live.rules.filter(r=>r.cap.startsWith('Color group')).map(r=>r.cap).sort()).toEqual(captions);
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
      expect(live.notice).toMatch(/Updated 3/);
      const notice=page.getByRole('status').filter({hasText:'specific-color caption'});
      await expect(notice).toHaveAttribute('aria-live','polite');
      await notice.scrollIntoViewIfNeeded();
      const noticeBounds=await notice.boundingBox();
      expect(noticeBounds.x).toBeGreaterThanOrEqual(0);
      expect(noticeBounds.x+noticeBounds.width).toBeLessThanOrEqual(width);
      await page.screenshot({path:testInfo.outputPath('caption-notice.png')});
      await history(page,'Undo');
      expect((await inspect(page)).rules).toEqual(before.rules);
      await history(page,'Redo');
      expect((await inspect(page)).rules).toEqual(live.rules);
      await generate(page);
      const generated=await inspect(page);
      expect(generated.result.filter(e=>e.caption.startsWith('Color group')).sort((a,b)=>a.caption.localeCompare(b.caption))).toEqual(live.result.filter(e=>e.caption.startsWith('Color group')).sort((a,b)=>a.caption.localeCompare(b.caption)));
      // Rename/recolor acts on the canonical source rows and leaves siblings intact.
      await page.evaluate(async()=>{
        const a=window.__GBDRAW_APP__;
        const index=a.manualSpecificRules.findIndex(r=>r.cap==='Color group [#112233]');
        await a.setSpecificRuleField(index,'cap','Primary complex');
        await a.updateLegendEntryColor(a.legendEntries.findIndex(e=>e.caption==='Primary complex'),'#224466');
        a.sortLegendEntries();
      });
      const edited=await inspect(page);
      expect(edited.rules.find(r=>r.val===source.ids[0]).cap).toBe('Primary complex');
      expect(edited.rules.find(r=>r.val===source.ids[0]).color).toBe('#224466');
      expect(edited.rules.find(r=>r.val===source.ids[1]).cap).toBe('Color group [#445566]');
      await generate(page);
      const final=await inspect(page);
      expect(final.result.find(e=>e.caption==='Primary complex').color).toBe('#224466');
      const pendingTsv=page.waitForEvent('download');
      await page.evaluate(()=>window.__GBDRAW_APP__.downloadSpecificRulesTsv());
      const tsvPath=testInfo.outputPath('canonical-rules.tsv');
      await (await pendingTsv).saveAs(tsvPath);
      const tsv=await fs.readFile(tsvPath,'utf8');
      expect(tsv).toContain('\tPrimary complex\n');
      expect(tsv).toContain('\tColor group [#445566]\n');
      expect(await page.evaluate(async()=>{const {readFileText}=await import('./js/services/file-content-cache.js');return readFileText(window.__GBDRAW_APP__.files.t_color);})).toBe(source.text);
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
 root=E.fromstring(result.drawing.tostring())
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
      await page.evaluate(async()=>{const a=window.__GBDRAW_APP__;await a.removeSpecificRule(a.manualSpecificRules.findIndex(r=>r.cap==='Color group [#445566]'));});
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


test('generated-caption collisions roll back and old Session drafts normalize on the next edit and Generate',async({browser},testInfo)=>{
  test.setTimeout(1_800_000);
  const page=await load(browser,'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
  page.setDefaultTimeout(180000);
  let fresh;
  try {
    await generate(page);
    await page.evaluate(async()=>{
      const a=window.__GBDRAW_APP__;
      a.newLegendCaption='Conflict [#112233]';a.newLegendColor='#abcdef';
      await a.addNewLegendEntry();
    });
    await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.legendEntries.some(e=>e.caption==='Conflict [#112233]'))).toBe(true);
    const ids=await page.evaluate(async()=>{
      const a=window.__GBDRAW_APP__,{getFeatureGenerationHash}=await import('./js/app/feature-utils.js');
      const ids=a.extractedFeatures.filter(f=>f.type==='CDS').slice(0,2).map(getFeatureGenerationHash);
      Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:ids[0],color:'#445566',cap:'Conflict'});
      await a.addSpecificRule();return ids;
    });
    const before=await inspect(page);
    const warningBefore=await page.evaluate(async()=>JSON.stringify((await import('./js/state.js')).state.annotationWarnings.value));
    await page.evaluate(async id=>{
      const a=window.__GBDRAW_APP__;
      Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:id,color:'#112233',cap:'Conflict'});
      await a.addSpecificRule();
    },ids[1]);
    const rejected=await inspect(page);
    expect(rejected.rules).toEqual(before.rules);expect(rejected.svg).toBe(before.svg);
    expect(rejected.result).toEqual(before.result);
    expect(await page.evaluate(async()=>JSON.stringify((await import('./js/state.js')).state.annotationWarnings.value))).toBe(warningBefore);
    // A historical draft fixture: ordinary saved strings with an unchanged preview.
    await page.evaluate(async ids=>{
      const a=window.__GBDRAW_APP__;
      const {manualSpecificRules,featureColorOverrides}=a;
      manualSpecificRules.splice(0,manualSpecificRules.length,
        ...ids.map((id,i)=>({feat:'CDS',qual:'hash',val:id,color:['#112233','#445566'][i],cap:'Historical'})));
      await a.renameLegendEntry(a.legendEntries.findIndex(entry=>entry.caption==='Conflict'),'Historical');
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
    const normalized=await inspect(fresh);
    expect(normalized.rules.map(r=>r.cap)).toEqual(['Historical [#112233]','Historical [#445566]']);
    expect(normalized.notice).toMatch(/Updated 2/);
    for(const [caption,color] of [['Historical [#112233]','#112233'],['Historical [#445566]','#445566']])
      expect(normalized.result.find(e=>e.caption===caption)?.color).toBe(color);
    expect(normalized.result.some(e=>e.caption==='Historical')).toBe(false);
    expect(normalized.result.find(e=>e.caption==='Conflict [#112233]').color).toBe('#abcdef');
    await generate(fresh);
    expect((await inspect(fresh)).result.filter(e=>e.caption.startsWith('Historical'))).toEqual(normalized.result.filter(e=>e.caption.startsWith('Historical')));
  } finally {await page.context().close();if(fresh)await fresh.context().close();}
});


test('Linear comparison keeps both legend orientations and their swatches canonical',async({browser},testInfo)=>{
  test.setTimeout(1_800_000);
  const page=await load(browser,'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
  page.setDefaultTimeout(180000);
  try {
    await generate(page);
    await page.evaluate(async()=>{
      const a=window.__GBDRAW_APP__,{getFeatureGenerationHash}=await import('./js/app/feature-utils.js');
      const features=a.extractedFeatures.filter(f=>f.type==='CDS').slice(0,2);
      for(let i=0;i<features.length;i++) {
        Object.assign(a.newSpecRule,{feat:'CDS',qual:'hash',val:getFeatureGenerationHash(features[i]),color:['#112233','#445566'][i],cap:'Comparison group'});
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
    await generate(page);
    const fresh=await inspect(page);
    expect(fresh.result.filter(e=>e.caption.startsWith('Comparison group'))).toEqual(live.result.filter(e=>e.caption.startsWith('Comparison group')));
    const pending=page.waitForEvent('download');
    await page.evaluate(()=>window.__GBDRAW_APP__.downloadSVG());
    await (await pending).saveAs(testInfo.outputPath('comparison.svg'));
    await page.locator('.origin-top svg').screenshot({path:testInfo.outputPath('comparison.svg.png')});
    await fs.writeFile(testInfo.outputPath('dual-legend.json'),JSON.stringify({...fresh,svg:undefined},null,2));
  }finally{await page.context().close();}
});
