"""Characterize S00 boundaries using unchanged runtime and disposable sessions."""
import argparse
import copy
import functools
import gzip
import hashlib
import json
import math
import platform
import subprocess
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from importlib.metadata import version
from pathlib import Path

from playwright.sync_api import sync_playwright
from gbdraw.session import load_session_document

parser = argparse.ArgumentParser()
parser.add_argument('--source-root', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--variants-only', action='store_true')
parser.add_argument('--only-case', action='append')
parser.add_argument('--only-variant', action='append')
args = parser.parse_args()
root = args.source_root.resolve()
repo = Path(__file__).resolve().parents[4]
cases = [] if args.variants_only else json.loads(Path(__file__).with_name('scalar-fixtures.json').read_text())
if args.only_case:
    cases = [sample for sample in cases if sample['name'] in args.only_case]
args.artifacts.mkdir(parents=True, exist_ok=True)
fixture = 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'
seed = json.loads((root / fixture).read_text())
source_sha = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
# A snapshot must match the handoff runtime, not a modified candidate.
source_paths = subprocess.check_output(['git', '-C', str(repo), 'ls-files', 'gbdraw/web/js', 'gbdraw/web/index.html'], text=True).splitlines()
for path in source_paths:
    assert (root / path).read_bytes() == (repo / path).read_bytes(), path


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()


def encode(value):
    if isinstance(value, float) and not math.isfinite(value):
        return {"$number": "NaN" if math.isnan(value) else "Infinity" if value > 0 else "-Infinity"}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value


def native(path):
    try:
        value = load_session_document(path).to_dict()
        row = next((s for s in value.get('config', {}).get('adv', {}).get('circular_track_slots', []) if s['id'] == 'gc_content'), None)
        if row is None:
            return {'accepted': True, 'configPresent': bool(value.get('config'))}
        return encode({'accepted': True, 'width': row['width'], 'radius': row['radius']})
    except Exception as error:
        return {'accepted': False, 'error': str(error)}


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


server = ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(QuietHandler, directory=str(root)))
threading.Thread(target=server.serve_forever, daemon=True).start()
url = f'http://127.0.0.1:{server.server_port}/gbdraw/web/index.html'
observations = []
variants = []
dialogs = []
page_errors = []
console_errors = []
external = []
read = """async () => {
 const s = (await import('./js/state.js')).state;
 const c = await import('./js/services/config.js');
 const row = s.adv.circular_track_slots.find(r => r.id === 'gc_content');
 const h = window.__GBDRAW_HISTORY__;
 const encode = v => JSON.parse(JSON.stringify(v, (_,x) => typeof x === 'number' && !Number.isFinite(x) ? {$number:String(x)} : x));
 return {width:row ? encode(row.width) : {$absent:true},radius:row ? encode(row.radius) : {$absent:true},request:c.getCommittedCanonicalRenderRequest(),result:s.results.value.map(r=>r.content),undo:h.getUndoCount(),redo:h.getRedoCount(),workers:window.__s00Workers.length};
}"""


def open_page(browser):
    page = browser.new_page(viewport={'width': 1440, 'height': 1000})
    page.set_default_timeout(180_000)
    page.on('dialog', lambda dialog: (dialogs.append(dialog.message), dialog.accept('S00 settings' if dialog.type == 'prompt' else None)))
    page.on('pageerror', lambda error: page_errors.append(str(error)))
    page.on('console', lambda message: console_errors.append(message.text) if message.type == 'error' else None)
    page.route('**/*', lambda route: route.continue_() if route.request.url.startswith(('http://127.0.0.1:', 'blob:', 'data:')) else (external.append(route.request.url), route.abort()))
    page.add_init_script("""window.__s00Workers=[]; const W=window.Worker; window.Worker=new Proxy(W,{construct(t,a){window.__s00Workers.push(String(a[0]));return Reflect.construct(t,a);}}); window.__s00Events=[]; window.__GBDRAW_TEST_HOOKS__={onSessionLifecycleEvent:e=>window.__s00Events.push(e)};""")
    page.goto(url)
    page.wait_for_function('window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions).length > 0')
    return page


def load(page, path):
    count = page.evaluate("window.__s00Events.filter(e=>e.name==='interactiveReady').length")
    start = len(dialogs)
    page.locator('input[type=file][accept^=".json,"]').set_input_files(path)
    page.wait_for_function("n => window.__s00Events.filter(e=>e.name==='interactiveReady').length>n && !window.__GBDRAW_APP__.sessionImportPending", arg=count)
    return {'accepted': not any(message.startswith('Failed to load session:') for message in dialogs[start:]), 'dialogs': dialogs[start:]}


try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        browser_version = browser.version
        page = reader = None
        if cases:
            page = open_page(browser)
            assert load(page, root / fixture)['accepted']
            baseline = page.evaluate(read)
            reader = open_page(browser)
            assert load(reader, root / fixture)['accepted']
        for sample in cases:
            name = sample['name']
            print('CHECK', name, flush=True)
            before = page.evaluate(read)
            history = page.evaluate("""async input => {
              const s = (await import('./js/state.js')).state;
              const h = window.__GBDRAW_HISTORY__;
              const decode = v => v && typeof v==='object' ? ('$number' in v ? Number(v.$number) : Object.fromEntries(Object.entries(v).map(([k,x])=>[k,decode(x)]))) : v;
              await h.runUndoable('S00 scalar fixture', async () => {
                const row=s.adv.circular_track_slots.find(r=>r.id==='gc_content');
                row.width=decode(input);row.radius=decode(input);await window.Vue.nextTick();
              });
              return {width:h.getCurrentIntent().config.adv.circular_track_slots.find(r=>r.id==='gc_content').width,radius:h.getCurrentIntent().config.adv.circular_track_slots.find(r=>r.id==='gc_content').radius};
            }""", sample['input'])
            after = page.evaluate(read)
            assert after['request'] == before['request'] and after['result'] == before['result']
            await_undo = page.evaluate('window.__GBDRAW_HISTORY__.undo()')
            undone = page.evaluate(read)
            await_redo = page.evaluate('window.__GBDRAW_HISTORY__.redo()')
            redone = page.evaluate(read)
            entry = {'name': name, 'input': sample['input'], 'draftAfterEdit': {k: after[k] for k in ('width', 'radius')},
                     'history': {'captured': history, 'undoReturned': await_undo, 'redoReturned': await_redo,
                                 'undoRestored': all(undone[k] == before[k] for k in ('width', 'radius')),
                                 'redoRestored': all(redone[k] == after[k] for k in ('width', 'radius')),
                                 'redone': {k: redone[k] for k in ('width', 'radius')}},
                     'requestAndResultPreserved': True, 'workersDuringEdit': after['workers'] - before['workers']}
            if not entry['history']['redoRestored']:
                # Probe writer admission of the ORIGINAL raw value, independently of History loss.
                page.evaluate("""async input => {const s=(await import('./js/state.js')).state; const decode=v=>v&&typeof v==='object'?('$number' in v?Number(v.$number):Object.fromEntries(Object.entries(v).map(([k,x])=>[k,decode(x)]))):v; const r=s.adv.circular_track_slots.find(r=>r.id==='gc_content');r.width=decode(input);r.radius=decode(input);await window.Vue.nextTick();}""", sample['input'])
            entry['draftAtSave'] = {key: value for key, value in page.evaluate(read).items() if key in ('width', 'radius')}
            if sample['valid']:
                assert entry['history']['undoRestored'] and entry['history']['redoRestored'], name
                with page.expect_download() as pending:
                    page.get_by_role('button', name='Save Session', exact=True).click()
                saved_path = args.artifacts / f'{name}.gbdraw-session.json.gz'
                pending.value.save_as(saved_path)
                saved = json.loads(gzip.decompress(saved_path.read_bytes()))
                saved_row = next(r for r in saved['config']['adv']['circular_track_slots'] if r['id'] == 'gc_content')
                assert saved_row['width'] == sample['input'] and saved_row['radius'] == sample['input'], name
                assert saved['renderRequest'] == baseline['request'], name
                entry['save'] = {'accepted': True, 'via': 'Save Session button and download', 'width': saved_row['width'], 'radius': saved_row['radius'], 'committedRequestPreserved': True, 'downloadSha256': hashlib.sha256(saved_path.read_bytes()).hexdigest()}
                entry['nativeSaved'] = native(saved_path)
                assert entry['nativeSaved']['accepted'], (name, entry['nativeSaved'])
                fresh = open_page(browser)
                entry['loadSaved'] = load(fresh, saved_path)
                restored = fresh.evaluate(read)
                entry['loadSaved'].update({'width': restored['width'], 'radius': restored['radius'], 'requestPreserved': restored['request'] == baseline['request'], 'resultPreserved': restored['result'] == after['result']})
                assert entry['loadSaved']['accepted'] and entry['loadSaved']['resultPreserved'] and entry['loadSaved']['requestPreserved'], name
                assert restored['width'] == sample['input'] and restored['radius'] == sample['input'], name
                fresh.close()
            else:
                # Execute the real writer; errors are measurements, never replaced with null.
                entry['save'] = page.evaluate("""async () => {try {const r=await (await import('./js/services/config.js')).exportSession('S00 invalid fixture');return {accepted:r.status==='saved',status:r.status};}catch(e){return {accepted:false,error:e.message};}}""")
            page.evaluate('window.__GBDRAW_HISTORY__.undo()')
            # Native/browser readers are independently probed against bytes, including malformed drafts.
            synthetic = copy.deepcopy(seed)
            row = next(r for r in synthetic['config']['adv']['circular_track_slots'] if r['id'] == 'gc_content')
            def decode_py(value):
                if isinstance(value, dict):
                    if '$number' in value:
                        return float(value['$number'])
                    return {key: decode_py(item) for key, item in value.items()}
                return value
            row['width'] = row['radius'] = decode_py(sample['input'])
            synthetic_path = args.artifacts / f'{name}.input.json'
            synthetic_path.write_text(json.dumps(synthetic))
            entry['nativeSynthetic'] = native(synthetic_path)
            reader_before = reader.evaluate(read)
            entry['loadSynthetic'] = load(reader, synthetic_path)
            reader_after = reader.evaluate(read)
            entry['loadSynthetic'].update({'width': reader_after['width'], 'radius': reader_after['radius'], 'requestPreserved': reader_after['request'] == reader_before['request'], 'resultPreserved': reader_after['result'] == reader_before['result']})
            if not entry['loadSynthetic']['accepted']:
                assert all(reader_after[k] == reader_before[k] for k in ('width','radius','request','result')), name
            observations.append(entry)
            (args.artifacts / 'matrix-progress.json').write_text(json.dumps(observations, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
        if args.variants_only:
            for variant in args.only_variant or ('disabled', 'inactive', 'settings-only', 'settings-only-control', 'cli-origin'):
                print('VARIANT', variant, flush=True)
                target = open_page(browser)
                if variant == 'cli-origin':
                    cli_path = args.artifacts / 'cli-origin.gbdraw-session.json'
                    native_cli = native(cli_path)
                    assert native_cli['accepted'] and not native_cli['configPresent']
                    admitted = load(target, cli_path)
                    assert admitted['accepted'], admitted
                elif variant not in ('settings-only','settings-only-control'):
                    assert load(target, root / fixture)['accepted']
                if variant not in ('cli-origin','settings-only-control'):
                    target.evaluate("""async ({variant,slots}) => {
                      const s=(await import('./js/state.js')).state;
                      if(variant==='settings-only') s.adv.circular_track_slots=slots.filter(r=>['features','gc_content'].includes(r.id));
                      const row=s.adv.circular_track_slots.find(r=>r.id==='gc_content');
                      row.width={value:'1.',unit:'px'}; row.radius={value:'1e-3',unit:'factor'};
                      if(variant==='disabled') row.enabled=false;
                      if(variant==='inactive') s.adv.circular_track_slots_enabled=false;
                      await window.Vue.nextTick();
                    }""", {'variant': variant, 'slots': seed['config']['adv']['circular_track_slots']})
                state_before = target.evaluate(read)
                downloads = []
                target.on('download', lambda download: downloads.append(download))
                count = target.evaluate("window.__s00Events.filter(e=>e.name==='session-save-pending-cleared').length")
                target.get_by_role('button', name='Save Session', exact=True).click()
                target.wait_for_function("n=>window.__s00Events.filter(e=>e.name==='session-save-pending-cleared').length>n", arg=count)
                if not downloads:
                    error = target.evaluate("window.__GBDRAW_APP__.errorLog")
                    print('SAVE REJECTED', variant, error, flush=True)
                    variants.append({'name':variant,'saveAccepted':False,'error':error,'consoleErrors':console_errors[:],'width':state_before['width'],'radius':state_before['radius'],'requestAndResultPreserved':all(target.evaluate(read)[k]==state_before[k] for k in ('request','result'))})
                    (args.artifacts / 'variant-progress.json').write_text(json.dumps(variants, indent=2, ensure_ascii=False) + '\n')
                    target.close()
                    continue
                path = args.artifacts / f'variant-{variant}.gbdraw-session.json.gz'
                downloads[0].save_as(path)
                saved = json.loads(gzip.decompress(path.read_bytes()))
                row = next((r for r in saved['config']['adv']['circular_track_slots'] if r['id']=='gc_content'), None)
                native_result = native(path)
                assert native_result['accepted'], native_result
                fresh = open_page(browser)
                loaded = load(fresh, path)
                assert loaded['accepted'], loaded
                restored = fresh.evaluate(read)
                assert all(restored[k]==state_before[k] for k in ('width','radius','request','result')), variant
                flags = fresh.evaluate("""() => ({enabled:window.__GBDRAW_APP__.adv.circular_track_slots.find(r=>r.id==='gc_content')?.enabled ?? null,useSlots:window.__GBDRAW_APP__.adv.circular_track_slots_enabled})""")
                if variant == 'disabled':
                    assert flags['enabled'] is False
                if variant == 'inactive':
                    assert flags['useSlots'] is False
                if variant == 'settings-only':
                    assert saved['renderRequest'] is None and not saved['results']
                variants.append({'name':variant,'saveAccepted':True,'loadAccepted':True,'native':native_result,'width':row['width'] if row else {'$absent':True},'radius':row['radius'] if row else {'$absent':True},'restoredFlags':flags,'requestAndResultPreserved':True,'sourceConfigPresent':False if variant=='cli-origin' else True,'downloadSha256':hashlib.sha256(path.read_bytes()).hexdigest()})
                (args.artifacts / 'variant-progress.json').write_text(json.dumps(variants, indent=2, ensure_ascii=False) + '\n')
                fresh.close()
                target.close()
        if reader is not None:
            reader.close()
        if page is not None:
            page.close()
        browser.close()
finally:
    server.shutdown()
    server.server_close()
result = {'sourceSha': source_sha, 'snapshotRuntimeMatchesCheckout': True,
          'fixture': fixture, 'fixtureSha256': hashlib.sha256((root / fixture).read_bytes()).hexdigest(),
          'runtimeTreeSha256': digest('\n'.join(f'{path}:{hashlib.sha256((root/path).read_bytes()).hexdigest()}' for path in source_paths)),
          'wheelSha256': hashlib.sha256((root / 'gbdraw/web/gbdraw-0.14.0-py3-none-any.whl').read_bytes()).hexdigest(),
          'pythonVersion': platform.python_version(), 'playwrightVersion': version('playwright'), 'browserVersion': browser_version,
          'observations': observations, 'variants': variants, 'pageErrors': page_errors, 'consoleErrors': console_errors, 'externalRequests': external,
          'limits': ['No new controls, selector/IME transaction, Generate geometry, SVG download/replay or Auto unit lifetime was implemented or tested.',
                     'Fixture scalar changes call the existing History manager in a real app; they do not prove future DOM input transaction wiring.',
                     'Native acceptance alone is not Web draft admission. NaN/Infinity JSON bytes are deliberately nonstandard negative reader inputs.']}
assert not page_errors, page_errors
assert not external, external
args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
print('WROTE', args.output, flush=True)
