import argparse,hashlib,json,pathlib,platform,subprocess
from playwright.sync_api import sync_playwright
parser=argparse.ArgumentParser(description='Observe issue 597 base behavior without calling Generate or refresh actions.')
parser.add_argument('--repo',type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[4])
parser.add_argument('--url',default='http://127.0.0.1:4597')
parser.add_argument('--output',type=pathlib.Path,default=pathlib.Path('/tmp/issue597-browser-evidence.json'))
args=parser.parse_args()
root=args.repo.resolve()
source=(root/'tests/fixtures/sessions/cli-web-mito.gb').read_bytes()
results={"baseSha":subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),"fixtureSha256":hashlib.sha256(source).hexdigest(),"fixture":"tests/fixtures/sessions/cli-web-mito.gb","platform":platform.platform()}
with sync_playwright() as p:
    browser=p.chromium.launch()
    results['browser']=browser.version
    page=browser.new_page(viewport={"width":1440,"height":1000})
    errors=[]
    page.on('pageerror',lambda e: errors.append(str(e)))
    page.on('dialog',lambda d:d.accept())
    page.add_init_script("window.__probeWorkers=[]; const W=window.Worker; window.Worker=class extends W { constructor(...a){super(...a); window.__probeWorkers.push(String(a[0]));} };")
    page.goto(args.url.rstrip('/')+'/gbdraw/web/index.html')
    page.wait_for_function("window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length > 0")
    upload=page.get_by_label('GenBank/DDBJ File',exact=True)
    results['multipleUploadEnabled']=upload.evaluate('(e)=>e.multiple')
    upload.set_input_files({"name":"native-multi.gb","mimeType":"text/plain","buffer":source+b'\n'+source})
    page.wait_for_function("window.__GBDRAW_APP__.circularRecordList.length===2")
    results['nativeBeforeGenerate']=page.evaluate("""async ()=>{const {state}=await import('./js/state.js');return {discoveryStatus:state.circularRecordDiscovery.status,records:state.circularRecordList.value,workers:window.__probeWorkers,resultCount:state.results.value.length,sharedCanvas:state.form.multi_record_canvas,cropPanelOpen:document.querySelector('[data-circular-record-presentation]').open};}""")
    results['nativeBeforeGenerate']['rotationSpinbuttons']=page.get_by_role('spinbutton',name='Display start',exact=False).count()
    results['nativeBeforeGenerate']['manualLoadLinks']=page.get_by_role('button',name='Load record rotation controls',exact=True).count()
    page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]').set_input_files(str(root/'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'))
    page.wait_for_function("window.__GBDRAW_APP__.sessionImportPending===false && window.__GBDRAW_APP__.results.length===1")
    results['savedPreviewAfterLoad']=page.evaluate("""async ()=>{const {state}=await import('./js/state.js');return {discoveryStatus:state.circularRecordDiscovery.status,recordCount:state.circularRecordList.value.length,workers:window.__probeWorkers,resultCount:state.results.value.length};}""")
    upload=page.get_by_label('GenBank/DDBJ File',exact=True)
    upload.set_input_files({"name":"replacement.gb","mimeType":"text/plain","buffer":source})
    page.wait_for_function("window.__GBDRAW_APP__.circularRecordList.length===1")
    results['replacementAfterLoad']=page.evaluate("""async ()=>{const {state}=await import('./js/state.js');return {discoveryStatus:state.circularRecordDiscovery.status,recordCount:state.circularRecordList.value.length,workers:window.__probeWorkers,resultCount:state.results.value.length};}""")
    results['pageErrors']=errors
    browser.close()
args.output.write_text(json.dumps(results,ensure_ascii=False,indent=2))
print(json.dumps(results,ensure_ascii=False,indent=2))
