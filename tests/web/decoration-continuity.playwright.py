"""Real Gallery/pointer acceptance for Issue 599 S01 (Python Playwright).

Run: python tests/web/decoration-continuity.playwright.py --evidence-dir /tmp/issue599-browser
"""
import argparse
import base64
import functools
import hashlib
import http.server
import importlib.metadata
import json
from pathlib import Path
import threading
import zipfile

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
READ = """async () => {
  const a=window.__GBDRAW_APP__, c=await import('./js/app/legend-layout/composition-actions.js'),
    {state}=await import('./js/state.js');
  const svg=a.svgContainer.querySelector('svg'), b=c.bindCompositionMetadata(svg);
  return {deltas:c.compositionUserDeltas(svg), primaryIds:b.primary.targets.map(t=>t.id),
    automatic:b.metadata.legend?.automaticTranslation,
    scaleText:svg.querySelector('#length_bar')?.textContent,
    biology:state.featureCatalog.value.items.map(i=>i.biologicalFeatures.map(f=>[f.recordKey,f.biologicalFeatureId,f.start,f.end,f.strand])),
    history:a.historyEntries?.length};
}"""


def assert_deltas(actual, expected):
    assert actual['primaryIds'] == expected['primaryIds']
    assert actual['biology'] == expected['biology']
    for role in ['legend', 'title', 'primary']:
        assert actual['deltas'][role] == expected['deltas'][role], (role, actual, expected)


def load_session(page, path):
    result=page.evaluate("""async ({name,encoded}) => {
      const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));
      const file=new File([bytes],name,{type:'application/octet-stream'});
      return await window.__GBDRAW_APP__.importSession({target:{files:[file],value:''}});
    }""", {'name':path.name,'encoded':base64.b64encode(path.read_bytes()).decode()})
    assert result['status']=='ok', result


def drag(page, selector, dx=24, dy=12):
    target = page.locator(f'[aria-label="Result Preview"] {selector}').first
    box = target.bounding_box()
    assert box, selector
    page.evaluate("""({x,y}) => {
      const a=window.__GBDRAW_APP__, v=document.querySelector('.preview-viewport').getBoundingClientRect();
      a.canvasPan.x+=v.x+v.width/2-x; a.canvasPan.y+=v.y+v.height/2-y;
    }""", {'x': box['x']+box['width']/2, 'y': box['y']+box['height']/2})
    page.wait_for_timeout(350)
    box = target.bounding_box()
    x, y = box['x']+box['width']/2, box['y']+box['height']/2
    hit = page.evaluate('({x,y,selector}) => Boolean(document.elementFromPoint(x,y)?.closest(selector.split(" ")[0]))', {'x':x, 'y':y, 'selector':selector})
    assert hit, selector
    page.mouse.move(x,y); page.mouse.down(); page.mouse.move(x+dx,y+dy,steps=6); page.mouse.up()
    page.wait_for_timeout(300)
    return {'selector': selector, 'intended': hit}


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *_args):
        pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence-dir', type=Path, required=True)
    args = parser.parse_args()
    args.evidence_dir.mkdir(parents=True, exist_ok=True)
    wheel = next((ROOT / 'gbdraw/web').glob('gbdraw-*.whl'))
    with zipfile.ZipFile(wheel) as archive:
        mismatches = [name for name in archive.namelist() if name.startswith('gbdraw/')
            and name.endswith(('.py', '.toml')) and (ROOT / name).exists()
            and archive.read(name) != (ROOT / name).read_bytes()]
    assert not mismatches, mismatches
    evidence = {'wheelSha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
        'wheelSourceMismatches': mismatches, 'playwrightVersion': importlib.metadata.version('playwright'),
        'generations': []}
    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0),
        functools.partial(QuietHandler, directory=str(ROOT)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            evidence['chromiumVersion'] = browser.version
            page = browser.new_page(viewport={'width': 1440, 'height': 1000}, accept_downloads=True)
            page.set_default_timeout(180_000)
            page.on('dialog', lambda dialog: dialog.accept())
            page.on('framenavigated', lambda frame: print('navigation:',frame.url,flush=True) if frame==page.main_frame else None)
            page.add_init_script("""window.__s01Metrics={}; window.__GBDRAW_TEST_HOOKS__={
              onStructuralMetric: m => { window.__s01Metrics[m.name]=(window.__s01Metrics[m.name]||0)+m.value; }
            };""")
            page.goto(f'http://127.0.0.1:{server.server_port}/gbdraw/web/index.html')
            page.wait_for_function('() => window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length')
            for mode, fixture in [('circular', 'HmmtDNA_basic_circular'), ('linear', 'lambda_basic_linear')]:
                load_session(page, ROOT / f'gbdraw/web/gallery/sessions/{fixture}.gbdraw-session.json')
                page.evaluate('() => { window.__s01Metrics={}; }')
                initial = page.evaluate("""async () => {
                  const a=window.__GBDRAW_APP__; a.form.plot_title='Composition continuity';
                  a.adv.plot_title_position='top'; return await a.runAnalysis();
                }""")
                assert initial['status'] == 'ok', initial
                zero_metrics=page.evaluate('() => window.__s01Metrics')
                assert zero_metrics.get('applicationSvgParseCount',0)==0, zero_metrics
                page.evaluate("""() => { const a=window.__GBDRAW_APP__; a.layoutRepositionMode=true; a.zoom=0.6; }""")
                selectors = ['#legend text', '#plot_title text'] + (['#length_bar text'] if mode == 'linear' else [])
                hits = []
                for selector in selectors:
                    hits.append(drag(page, selector))
                dragged = page.evaluate(READ)
                for role in ['legend', 'title']:
                    assert all(abs(v-e) < 1e-8 for v,e in zip(dragged['deltas'][role], [40,20])), dragged
                if mode == 'linear':
                    idx = dragged['primaryIds'].index('length_bar')
                    assert all(abs(v-e)<1e-8 for v,e in zip(dragged['deltas']['primary'][idx], [40,20])), dragged
                outcomes=[]
                for _ in range(2):
                    page.evaluate('() => { window.__s01Metrics={}; }')
                    outcome=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
                    assert outcome['status']=='ok', outcome
                    after=page.evaluate(READ); assert_deltas(after, dragged)
                    metrics=page.evaluate('() => window.__s01Metrics')
                    assert metrics.get('applicationSvgParseCount')==1, metrics
                    assert metrics.get('svgSerializationCount')==1, metrics
                    assert metrics.get('previewBinderInvocationCount')==1, metrics
                    assert metrics.get('canonicalCandidateExecutionCount')==1, metrics
                    assert metrics.get('generatedArtifactHeavyTraversalCount',0)==0, metrics
                    outcomes.append({'outcome':outcome, 'composition':after, 'metrics':metrics})
                # Color/font/title/side changes must use the new automatic baseline.
                changed=page.evaluate("""async () => {
                  const a=window.__GBDRAW_APP__; a.form.plot_title='A longer title after editing';
                  a.adv.plot_title_font_size=38; a.adv.legend_font_size=20;
                  a.adv.block_stroke_color='#223344'; a.form.legend='left';
                  return await a.runAnalysis();
                }""")
                assert changed['status']=='ok', changed
                changed_composition=page.evaluate(READ); assert_deltas(changed_composition, dragged)
                assert changed_composition['automatic'] != dragged['automatic']
                assert changed_composition['scaleText'] == dragged['scaleText']
                page.evaluate('async () => await window.__GBDRAW_APP__.undoHistory()')
                assert_deltas(page.evaluate(READ), dragged)
                page.evaluate('async () => await window.__GBDRAW_APP__.redoHistory()')
                assert_deltas(page.evaluate(READ), dragged)
                with page.expect_download() as download:
                    page.evaluate('async () => await window.__GBDRAW_APP__.downloadSVG()')
                export_path=args.evidence_dir/f'{mode}.svg'; download.value.save_as(export_path)
                exported=page.evaluate("""async text => {
                  const c=await import('./js/app/legend-layout/composition-actions.js');
                  return c.compositionUserDeltas(new DOMParser().parseFromString(text,'image/svg+xml').documentElement);
                }""", export_path.read_text())
                assert exported==dragged['deltas'], exported
                figure=browser.new_page(viewport={'width':1440,'height':1200})
                figure.set_content(export_path.read_text())
                figure.locator('svg').evaluate("svg => { svg.style.width='1300px'; svg.style.height='auto'; }")
                figure.locator('svg').screenshot(path=str(args.evidence_dir/f'{mode}-figure.png'))
                figure.close()
                with page.expect_download() as download:
                    page.evaluate("async () => { window.__GBDRAW_APP__.sessionTitle='S01 continuity'; await window.__GBDRAW_APP__.saveSessionWithTitle(); }")
                session_path=args.evidence_dir/f'{mode}.gbdraw-session.json'; download.value.save_as(session_path)
                load_session(page,session_path)
                assert_deltas(page.evaluate(READ), dragged)
                loaded=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
                assert loaded['status']=='ok', loaded
                assert_deltas(page.evaluate(READ), dragged)
                page.evaluate('async () => await window.__GBDRAW_APP__.resetLegendPosition()')
                reset=page.evaluate(READ)
                assert reset['deltas']['legend']==[0,0]
                assert reset['deltas']['title']==dragged['deltas']['title']
                page.evaluate('async () => await window.__GBDRAW_APP__.undoHistory()')
                assert_deltas(page.evaluate(READ),dragged)
                page.evaluate('async () => await window.__GBDRAW_APP__.redoHistory()')
                assert page.evaluate(READ)['deltas']['legend']==[0,0]
                page.evaluate('async () => await window.__GBDRAW_APP__.undoHistory()')
                assert_deltas(page.evaluate(READ),dragged)
                page.screenshot(path=str(args.evidence_dir/f'{mode}-decorated.png'))
                # Nonzero target removal is atomic and recoverable through Reset.
                page.evaluate("async () => { window.__GBDRAW_APP__.form.legend='none'; await window.Vue.nextTick(); }")
                page.wait_for_timeout(350)
                before=page.evaluate('() => window.__GBDRAW_APP__.results.map(r=>r.content)')
                failed=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
                assert failed['status']=='error', failed
                assert page.evaluate('() => window.__GBDRAW_APP__.results.map(r=>r.content)') == before
                assert 'Reset' in json.dumps(failed)
                page.evaluate('async () => await window.__GBDRAW_APP__.resetLayout()')
                recovered=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
                assert recovered['status']=='ok', recovered
                page.screenshot(path=str(args.evidence_dir/f'{mode}.png'))
                evidence['generations'].append({'mode':mode, 'fixture':fixture, 'pointerHits':hits,
                    'dragged':dragged, 'zeroMetrics':zero_metrics, 'repeat':outcomes, 'changed':changed_composition,
                    'sessionGenerate':loaded, 'targetFailure':failed, 'recovery':recovered})
                print(f'{mode}: pointer drag, Generate x2, baseline change, History, export, Session, failure/Reset PASS', flush=True)
            # Two real batch outputs; the second edit must survive while unselected.
            load_session(page,ROOT / 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')
            created=page.evaluate("""async () => {
              const a=window.__GBDRAW_APP__, {readFileText}=await import('./js/services/file-content-cache.js'),
                text=await readFileText(a.files.c_gb);
              a.files.c_gb=new File([text,'\\n',text], 'batch.gb', {type:'text/plain'});
              a.adv.circular_grouping_intent='batch'; a.form.multi_record_canvas=false;
              a.form.circular_record_selector=''; a.form.plot_title='Batch continuity';
              a.adv.plot_title_position='top'; await a.refreshCircularRecordOrder();
              return await a.runAnalysis();
            }""")
            assert created['status']=='ok', created
            assert page.evaluate('() => window.__GBDRAW_APP__.results.length')==2
            expected=[]
            for index, offset in [(0,(24,12)), (1,(-12,6))]:
                page.evaluate('async index => { const a=window.__GBDRAW_APP__; a.selectedResultIndex=index; await window.Vue.nextTick(); }',index)
                page.wait_for_timeout(500)
                page.evaluate('() => { const a=window.__GBDRAW_APP__; a.layoutRepositionMode=true; a.zoom=0.6; }')
                drag(page,'#legend text',*offset); drag(page,'#plot_title text',*offset)
                expected.append(page.evaluate(READ))
            page.evaluate("async () => { const a=window.__GBDRAW_APP__; a.selectedResultIndex=0; a.form.prefix='renamed-batch'; await window.Vue.nextTick(); }")
            page.wait_for_timeout(500)
            batch_runs=[]
            for _ in range(2):
                page.evaluate('() => { window.__s01Metrics={}; }')
                outcome=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
                assert outcome['status']=='ok', outcome
                metrics=page.evaluate('() => window.__s01Metrics')
                assert metrics.get('applicationSvgParseCount')==3, metrics # one old unselected + two candidates
                assert metrics.get('svgSerializationCount')==2, metrics
                assert metrics.get('previewBinderInvocationCount')==1, metrics
                assert metrics.get('canonicalCandidateExecutionCount')==1, metrics
                for index in range(2):
                    page.evaluate('async i => { window.__GBDRAW_APP__.selectedResultIndex=i; await window.Vue.nextTick(); }',index)
                    page.wait_for_timeout(500)
                    assert_deltas(page.evaluate(READ),expected[index])
                batch_runs.append(metrics)
            with page.expect_download() as download:
                page.evaluate("async () => { window.__GBDRAW_APP__.sessionTitle='S01 batch'; await window.__GBDRAW_APP__.saveSessionWithTitle(); }")
            batch_session=args.evidence_dir/'batch.gbdraw-session.json'; download.value.save_as(batch_session)
            load_session(page,batch_session)
            outcome=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
            assert outcome['status']=='ok', outcome
            for index in range(2):
                page.evaluate('async i => { window.__GBDRAW_APP__.selectedResultIndex=i; await window.Vue.nextTick(); }',index)
                page.wait_for_timeout(500); assert_deltas(page.evaluate(READ),expected[index])
            evidence['batch']={'expected':expected,'repeatMetrics':batch_runs,'sessionGenerate':outcome}
            print('batch: distinct pointer offsets, unselected Result, prefix, Generate x2, Session PASS',flush=True)
            browser.close()
    finally:
        server.shutdown(); server.server_close()
        (args.evidence_dir/'observations.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n')


if __name__ == '__main__':
    main()
