"""Issue 599 S02 real Gallery, gestures, restoration and export acceptance.

Run: python tests/web/layout-affordance.playwright.py --evidence-dir /tmp/issue599-s02-browser
"""
import argparse
import functools
import gzip
import hashlib
import http.server
import importlib.util
import json
import re
from pathlib import Path
import threading
import zipfile
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('continuity', Path(__file__).with_name('decoration-continuity.playwright.py'))
continuity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(continuity)
TOKEN = 'gbdraw-preview-layout-target'


def clean(text):
    assert TOKEN not in text
    assert 'Turn on Layout edit' not in text
    assert 'outline: 1px dashed #94a3b8' not in text
    # Standalone interactive CSS legitimately supplies its own cursors.
    assert 'cursor: help' not in text
    assert not re.search(r'<g\b[^>]*style=["\'][^"\']*cursor\s*:', text)


def point(page, selector):
    target = page.locator(f'.gbdraw-preview-surface {selector}').first
    box = target.bounding_box()
    assert box, selector
    page.evaluate("""({x,y}) => {
      const a=window.__GBDRAW_APP__, v=a.canvasContainerRef.getBoundingClientRect();
      a.canvasPan.x+=v.x+v.width/2-x; a.canvasPan.y+=v.y+v.height/2-y;
    }""", {'x': box['x']+box['width']/2, 'y': box['y']+box['height']/2})
    page.wait_for_timeout(350)
    box = target.bounding_box()
    p = {'x':box['x']+box['width']/2, 'y':box['y']+box['height']/2}
    assert page.evaluate('({p,s}) => Boolean(document.elementFromPoint(p.x,p.y)?.closest(s))', {'p':p, 's':selector}), selector
    return p


def status(page):
    return page.evaluate("""async () => {
      const a=window.__GBDRAW_APP__, c=await import('./js/app/legend-layout/composition-actions.js');
      return {deltas:c.compositionUserDeltas(a.svgContainer.querySelector('svg')),
        pan:{...a.canvasPan}, mode:a.layoutRepositionMode,
        history:window.__GBDRAW_HISTORY__.getUndoCount(), metrics:{...window.__s02Metrics},
        results:a.results.map(r=>r.content)};
    }""")


def bound(page, count):
    assert page.locator(f'.gbdraw-preview-surface .{TOKEN}').count() == count
    assert page.get_by_role('button',name='Toggle layout edit mode').get_attribute('aria-describedby') == 'preview-layout-help'


def cursor(page, selector):
    return page.locator(f'.gbdraw-preview-surface {selector}').first.evaluate('(e) => getComputedStyle(e).cursor')


def move(page, p, selector=None):
    page.mouse.move(p['x'],p['y']); page.mouse.down()
    if selector:
        assert cursor(page,selector)=='grabbing', (selector,cursor(page,selector))
    page.mouse.move(p['x']+24,p['y']+12,steps=6); page.mouse.up()
    page.wait_for_timeout(350)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--evidence-dir',type=Path,required=True)
    args=parser.parse_args(); args.evidence_dir.mkdir(parents=True,exist_ok=True)
    wheel=next((ROOT/'gbdraw/web').glob('gbdraw-*.whl'))
    with zipfile.ZipFile(wheel) as z:
        mismatches=[n for n in z.namelist() if n.startswith('gbdraw/') and n.endswith(('.py','.toml')) and (ROOT/n).exists() and z.read(n)!=(ROOT/n).read_bytes()]
    assert not mismatches, mismatches
    evidence={'wheelSha256':hashlib.sha256(wheel.read_bytes()).hexdigest(),'wheelSourceMismatches':mismatches,'modes':[]}
    server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(continuity.QuietHandler,directory=str(ROOT)))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    try:
        with sync_playwright() as pw:
            browser=pw.chromium.launch(); evidence['chromiumVersion']=browser.version
            page=browser.new_page(viewport={'width':1440,'height':1000},accept_downloads=True)
            page.set_default_timeout(180000); page.on('dialog',lambda d:d.accept())
            page.add_init_script("""window.__s02Metrics={}; window.__GBDRAW_TEST_HOOKS__={onStructuralMetric:m=>{
              window.__s02Metrics[m.name]=(window.__s02Metrics[m.name]||0)+m.value;
            }};""")
            url=f'http://127.0.0.1:{server.server_port}/gbdraw/web/index.html'
            page.goto(url); page.wait_for_function('() => window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions||{}).length')
            for mode,fixture in [('circular','HmmtDNA_basic_circular'),('linear','lambda_basic_linear')]:
                continuity.load_session(page,ROOT/f'gbdraw/web/gallery/sessions/{fixture}.gbdraw-session.json')
                generated=page.evaluate("""async () => { const a=window.__GBDRAW_APP__;
                  a.form.plot_title='Preview layout editing'; a.adv.plot_title_position='top';
                  return await a.runAnalysis(); }"""); assert generated['status']=='ok',generated
                page.evaluate('() => { const a=window.__GBDRAW_APP__; a.layoutRepositionMode=false; a.zoom=0.6; }')
                selectors=['#legend text','#plot_title text']+(['#length_bar text'] if mode=='linear' else [])
                bound(page,len(selectors)); operations=[]
                for selector in selectors:
                    p=point(page,selector); page.mouse.move(p['x'],p['y'])
                    assert cursor(page,selector)=='help', (selector,cursor(page,selector))
                    assert page.locator('.preview-layout-hint').is_visible()
                    assert page.locator(f'.gbdraw-preview-surface {selector.split()[0]}').evaluate('(e) => getComputedStyle(e).outlineStyle')=='dashed'
                    before=status(page)
                    for _ in range(4):
                        page.mouse.move(p['x']+1,p['y']+1); page.mouse.move(p['x'],p['y'])
                    hover=status(page); assert hover==before, (hover,before)
                    page.screenshot(path=str(args.evidence_dir/f'{mode}-{selector.split()[0][1:]}-hover.png'))
                    move(page,p); after=status(page)
                    assert after['deltas']==before['deltas'] and not after['mode']
                    assert after['history']==before['history'] and after['results']==before['results']
                    assert abs(after['pan']['x']-before['pan']['x']-24)<1e-8, (mode,selector,before['pan'],after['pan'])
                    assert abs(after['pan']['y']-before['pan']['y']-12)<1e-8
                    operations.append({'selector':selector,'offPan':[24,12],'hoverNoMutation':True})
                toggle=page.get_by_role('button',name='Toggle layout edit mode')
                toggle.focus(); assert 'Turn on Layout edit' in page.locator('#preview-layout-help').inner_text()
                page.keyboard.press('Space'); assert toggle.get_attribute('aria-pressed')=='true'
                assert not page.locator('.preview-layout-hint').count()
                for selector in selectors:
                    p=point(page,selector); assert cursor(page,selector)=='grab'
                    before=status(page); move(page,p,selector); after=status(page)
                    assert after['history']==before['history']+1,(before['history'],after['history'])
                    assert after['pan']==before['pan'] and after['mode']
                    assert cursor(page,selector)=='grab'
                    page.evaluate('async () => await window.__GBDRAW_APP__.undoHistory()'); bound(page,len(selectors))
                    assert status(page)['deltas']==before['deltas']
                    page.evaluate('async () => await window.__GBDRAW_APP__.redoHistory()'); bound(page,len(selectors))
                    assert status(page)['deltas']==after['deltas']

                # Editable labels keep their text cursor and take precedence over layout drag.
                page.evaluate('() => window.__GBDRAW_APP__.syncLabelEditor()')
                label_selector='text[data-label-editable="true"]'
                assert page.locator('.gbdraw-preview-surface '+label_selector).count()>0
                p=point(page,label_selector); before_label=status(page)
                assert cursor(page,label_selector)=='text'
                move(page,p)
                after_label=status(page)
                assert after_label['deltas']==before_label['deltas']
                assert after_label['pan']==before_label['pan']
                assert after_label['history']==before_label['history']
                close_popup=page.get_by_role('button',name='Close feature popup',exact=True)
                if close_popup.is_visible(): close_popup.click()

                # Shift keeps selection ownership; Ctrl retains the existing decoration gesture.
                p=point(page,'#plot_title text'); before=status(page)
                page.keyboard.down('Shift'); move(page,p); page.keyboard.up('Shift')
                assert status(page)['deltas']==before['deltas']
                p=point(page,'#plot_title text'); before=status(page)
                page.keyboard.down('Control'); move(page,p,'#plot_title text'); page.keyboard.up('Control')
                assert status(page)['history']==before['history']+1
                page.evaluate('async () => await window.__GBDRAW_APP__.undoHistory()')
                # Keyboard Enter returns to OFF before regeneration/export.
                toggle.focus(); page.keyboard.press('Enter'); assert toggle.get_attribute('aria-pressed')=='false'
                page.evaluate('() => window.__s02Metrics={}')
                before=status(page)
                outcome=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()'); assert outcome['status']=='ok',outcome
                bound(page,len(selectors)); after=status(page); assert after['deltas']==before['deltas']
                assert after['metrics'].get('previewBinderInvocationCount')==1,after['metrics']
                assert after['metrics'].get('canonicalCandidateExecutionCount')==1,after['metrics']
                for content in after['results']: clean(content)
                # Observe actual source handed to SVG serialization, PNG Image and PDF staging.
                page.evaluate("""() => {
                  window.__s02ExportSources=[];
                  const original=XMLSerializer.prototype.serializeToString;
                  XMLSerializer.prototype.serializeToString=function(node) {
                    const text=original.call(this,node);
                    if(node.localName==='svg') window.__s02ExportSources.push(text);
                    return text;
                  };
                }""")
                exports=[]
                for method,extension in [('downloadSVG','svg'),('downloadInteractiveSVG','interactive.svg'),('downloadPNG','png'),('downloadPDF','pdf')]:
                    p=point(page,'#legend text'); page.mouse.move(p['x'],p['y']); assert page.locator('.preview-layout-hint').is_visible()
                    page.evaluate('() => window.__s02ExportSources=[]')
                    with page.expect_download() as pending:
                        outcome=page.evaluate('async method => await window.__GBDRAW_APP__[method]()',method)
                    assert pending.value.failure() is None
                    path=args.evidence_dir/f'{mode}.{extension}'; pending.value.save_as(path)
                    assert path.stat().st_size>100
                    sources=page.evaluate('() => window.__s02ExportSources'); assert sources,method
                    for source in sources: clean(source)
                    (args.evidence_dir/f'{mode}-{method}-source.svg').write_text(sources[0])
                    if extension.endswith('svg'): clean(path.read_text())
                    if extension=='png': assert path.read_bytes().startswith(b'\x89PNG')
                    if extension=='pdf': assert path.read_bytes().startswith(b'%PDF')
                    exports.append({'method':method,'bytes':path.stat().st_size,'cleanSources':len(sources),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
                with page.expect_download() as pending:
                    page.evaluate("async () => { const a=window.__GBDRAW_APP__; a.sessionTitle='S02 affordance'; await a.saveSessionWithTitle(); }")
                path=args.evidence_dir/f'{mode}.gbdraw-session.json'; pending.value.save_as(path)
                clean(gzip.decompress(path.read_bytes()).decode())
                continuity.load_session(page,path); bound(page,len(selectors))
                restored=status(page); assert restored['deltas']==after['deltas']
                page.evaluate('() => window.__s02Metrics={}')
                loaded_generate=page.evaluate('async () => await window.__GBDRAW_APP__.runAnalysis()')
                assert loaded_generate['status']=='ok',loaded_generate
                bound(page,len(selectors)); post_load=status(page)
                assert post_load['deltas']==restored['deltas']
                assert post_load['metrics'].get('previewBinderInvocationCount')==1,post_load['metrics']
                assert post_load['metrics'].get('canonicalCandidateExecutionCount')==1,post_load['metrics']

                # One post-load drag still creates exactly one History action.
                toggle=page.get_by_role('button',name='Toggle layout edit mode'); toggle.click()
                assert toggle.get_attribute('aria-pressed')=='true'
                p=point(page,'#legend text'); before=status(page); move(page,p,'#legend text')
                assert status(page)['history']==before['history']+1
                evidence['modes'].append({'mode':mode,'operations':operations,'generateMetrics':after['metrics'],'exports':exports,'restoredDeltas':restored['deltas'],'postLoadGenerateMetrics':post_load['metrics']})
                print(mode+': OFF hover/pan, keyboard ON/drag, modifiers, Undo/Redo, Generate, Session, SVG/interactive/PNG/PDF PASS',flush=True)
            # Touch has the same persistent explanation and native toggle, without hover.
            touch=browser.new_context(viewport={'width':390,'height':844},has_touch=True,is_mobile=True,accept_downloads=True)
            tp=touch.new_page(); tp.on('dialog',lambda d:d.accept()); tp.goto(url)
            tp.wait_for_function('() => window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions||{}).length')
            continuity.load_session(tp,ROOT/'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json')
            help_text=tp.locator('#preview-layout-help'); help_text.scroll_into_view_if_needed(); assert help_text.is_visible()
            assert help_text.evaluate("""e => {
              const r=e.getBoundingClientRect(), x=r.x+r.width/2, y=r.y+r.height/2;
              return x>=0 && x<innerWidth && y>=0 && y<innerHeight
                && Boolean(document.elementFromPoint(x,y)?.closest('#preview-layout-help'));
            }""")
            tp.screenshot(path=str(args.evidence_dir/'touch-help.png'))

            tb=tp.get_by_role('button',name='Toggle layout edit mode'); tb.scroll_into_view_if_needed(); tb.tap()
            assert tb.get_attribute('aria-pressed')=='true'; tb.tap(); assert tb.get_attribute('aria-pressed')=='false'
            tp.screenshot(path=str(args.evidence_dir/'touch.png'))
            evidence['touch']={'viewport':[390,844],'helpText':help_text.inner_text(),'tapOnOff':True}
            touch.close();browser.close()
    finally:
        server.shutdown();server.server_close()
        (args.evidence_dir/'observations.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__': main()
