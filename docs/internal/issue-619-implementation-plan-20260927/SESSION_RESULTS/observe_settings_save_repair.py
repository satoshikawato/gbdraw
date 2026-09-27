"""Verify the settings-only writer repair using real downloads and fresh imports."""
import argparse
import copy
import functools
import gzip
import hashlib
import json
import platform
import subprocess
import threading
import zipfile
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from importlib.metadata import version
from pathlib import Path

from playwright.sync_api import sync_playwright
from gbdraw.session import load_session_document

parser = argparse.ArgumentParser()
parser.add_argument('--source-root', type=Path, required=True)
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
repo = Path(__file__).resolve().parents[4]
root = args.source_root.resolve()
args.artifacts.mkdir(parents=True, exist_ok=True)
fixture = 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'
seed = json.loads((root / fixture).read_text())
samples = json.loads(Path(__file__).with_name('scalar-fixtures.json').read_text())
paths = subprocess.check_output(['git', '-C', str(repo), 'ls-files', 'gbdraw/web/js', 'gbdraw/web/index.html', 'gbdraw'], text=True).splitlines()
paths = sorted(set(paths))
for path in paths:
    assert (root / path).read_bytes() == (repo / path).read_bytes(), path


def sha(value):
    return hashlib.sha256(value).hexdigest()


wheel = root / 'gbdraw/web/gbdraw-0.14.0-py3-none-any.whl'
with zipfile.ZipFile(wheel) as archive:
    python_paths = [path for path in paths if path.endswith('.py')]
    for path in python_paths:
        assert archive.read(path) == (root / path).read_bytes(), path

read = """async () => {
 const s=(await import('./js/state.js')).state;
 const c=await import('./js/services/config.js');
 const h=window.__GBDRAW_HISTORY__;
 const digest=async text=>Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(text)))).map(x=>x.toString(16).padStart(2,'0')).join('');
 const snapshot={config:c.buildConfigData(),request:c.getCommittedCanonicalRenderRequest(),
  resultHashes:await Promise.all(s.results.value.map(r=>digest(r.content))),
  resultNames:s.results.value.map(r=>r.name),
  metadata:{annotationWarnings:s.annotationWarnings.value,trackSlotGeometry:s.trackSlotResolvedGeometry.value},
  biological:(await import('./js/services/session-authority.js')).hasBiologicalSessionInputs({...s.files,linearSeqs:s.linearSeqs}),
  history:{undo:h.getUndoCount(),redo:h.getRedoCount()},workers:window.__repairWorkers.length};
 return JSON.parse(JSON.stringify(snapshot));
}"""
errors, external, dialogs, observations = [], [], [], []


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


server = ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(QuietHandler, directory=str(root)))
thread = threading.Thread(target=server.serve_forever, daemon=True)
thread.start()
url = f'http://127.0.0.1:{server.server_port}/gbdraw/web/index.html'


def fresh(browser):
    page = browser.new_page(viewport={'width': 1440, 'height': 1000})
    page.set_default_timeout(180_000)
    page.on('dialog', lambda dialog: (dialogs.append(dialog.message), dialog.accept('' if dialog.type == 'prompt' else None)))
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.route('**/*', lambda route: route.continue_() if route.request.url.startswith(('http://127.0.0.1:', 'blob:', 'data:')) else (external.append(route.request.url), route.abort()))
    page.add_init_script("""window.__repairEvents=[];window.__repairWorkers=[];
     window.__GBDRAW_TEST_HOOKS__={onSessionLifecycleEvent:e=>window.__repairEvents.push(e)};
     const W=window.Worker;window.Worker=new Proxy(W,{construct(t,a){window.__repairWorkers.push(String(a[0]));return Reflect.construct(t,a);}});""")
    page.goto(url)
    page.wait_for_function('window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions).length>0')
    return page


def load(page, path, accepted=True):
    count = page.evaluate("window.__repairEvents.filter(e=>e.name==='interactiveReady').length")
    start = len(dialogs)
    page.locator('input[type=file][accept^=".json,"]').set_input_files(path)
    page.wait_for_function("n=>window.__repairEvents.filter(e=>e.name==='interactiveReady').length>n && !window.__GBDRAW_APP__.sessionImportPending", arg=count)
    messages = dialogs[start:]
    success = messages == ['Session loaded successfully!']
    assert success == accepted, messages
    return {'accepted': success, 'dialogs': messages}


def scalar_pair(snapshot):
    row = next(r for r in snapshot['config']['adv']['circular_track_slots'] if r['id'] == 'gc_content')
    return {field: row[field] for field in ('width', 'radius')}


def save_roundtrip(browser, page, name, settings, expected_scalar=None, canonical=None):
    before = page.evaluate(read)
    assert page.evaluate(read) == before  # Read-only observation is inert.
    downloads = []
    page.on('download', lambda download: downloads.append(download))
    count = page.evaluate("window.__repairEvents.filter(e=>e.name==='session-save-pending-cleared').length")
    page.get_by_role('button', name='Save Session', exact=True).click()
    page.wait_for_function("n=>window.__repairEvents.filter(e=>e.name==='session-save-pending-cleared').length>n", arg=count)
    assert len(downloads) == 1, page.evaluate('window.__GBDRAW_APP__.errorLog')
    path = args.artifacts / f'{name}.gbdraw-session.json.gz'
    downloads[0].save_as(path)
    document = json.loads(gzip.decompress(path.read_bytes()))
    assert document['version'] == 44
    assert document['config'] == before['config']
    after_save = page.evaluate(read)
    if after_save != before:
        (args.artifacts / f'{name}.save-state-diff.json').write_text(json.dumps({'before': before, 'after': after_save}, indent=2))
        print('SAVE STATE DIFF', {key: (before[key], after_save[key]) for key in before if before[key] != after_save[key]}, flush=True)
    assert after_save == before
    if settings:
        assert document['renderRequest'] is None
        assert document['results'] == [] and document['editorState']['featureCatalog'] is None
        assert document['runMetadata'] == {} and not before['biological']
    else:
        assert document['renderRequest'] == before['request']
        assert document['runMetadata'] == before['metadata']
        assert document['resources']  # Biological bytes are retained.
    native = load_session_document(path).to_dict()
    assert native['config'] == document['config']
    assert native['renderRequest'] == document['renderRequest']
    assert native['runMetadata'] == document['runMetadata']
    restored_page = fresh(browser)
    try:
        admission = load(restored_page, path)
        after = restored_page.evaluate(read)
        for key in ('request', 'resultHashes', 'resultNames', 'biological'):
            assert after[key] == before[key], (name, key)
        if name == 'inactive-biological':
            # This control tests biological classification across modes. Current
            # Load selects the committed request mode; record that continuation.
            assert after['config']['adv']['circular_track_slots'] == before['config']['adv']['circular_track_slots']
        else:
            assert after['config'] == before['config'], (name, 'config')
        admission.update(configEqual=after['config'] == before['config'],
                         savedModeProfile=before['config']['modeProfiles']['activeMode'],
                         restoredModeProfile=after['config']['modeProfiles']['activeMode'])
        assert after['workers'] == before['workers']
        if settings:
            assert after['workers'] == 0
        if not settings:
            assert after['metadata'] == before['metadata']
        if expected_scalar is not None:
            assert scalar_pair(before) == expected_scalar and scalar_pair(after) == expected_scalar
            row = next(r for r in native['config']['adv']['circular_track_slots'] if r['id'] == 'gc_content')
            assert {k: row[k] for k in ('width', 'radius')} == expected_scalar
            for target in (page, restored_page):
                projected = target.evaluate("""async () => {
                  const s=(await import('./js/state.js')).state;
                  const p=(await import('./js/app/circular-track-slots.js')).buildCircularTrackSlotPayload(s.adv.circular_track_slots.find(r=>r.id==='gc_content'));
                  return {width:p.width,radius:p.radius};
                }""")
                assert projected == {'width': canonical, 'radius': canonical}, (name, projected)
        entry = {'name': name, 'settingsOnly': settings, 'save': {'accepted': True, 'via': 'Save Session button', 'download': str(path), 'sha256': sha(path.read_bytes())},
                 'freshLoad': admission, 'nativeReader': {'accepted': True, 'configRequestMetadataEqual': True},
                 'readAndSaveStateUnchanged': True, 'resultHashes': before['resultHashes'],
                 'history': before['history'], 'metadata': document['runMetadata'], 'freshLoadWorkerCount': after['workers']}
        if expected_scalar is not None:
            entry.update(rawScalars=expected_scalar, canonicalScalars={'width': canonical, 'radius': canonical})
        observations.append(entry)
        (args.artifacts / 'accepted-progress.json').write_text(json.dumps(observations, indent=2, ensure_ascii=False) + '\n')
        return document
    finally:
        restored_page.close()


try:
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        browser_version = browser.version
        page = fresh(browser)
        save_roundtrip(browser, page, 'fresh-control', True)
        page.close()
        for sample in (s for s in samples if s['valid']):
            print('SETTINGS', sample['name'], flush=True)
            page = fresh(browser)
            page.evaluate("""async ({slots,value})=>{
             const s=(await import('./js/state.js')).state;
             await window.__GBDRAW_HISTORY__.runUndoable('Settings scalar fixture',async()=>{
             s.adv.circular_track_slots=slots.filter(r=>['features','gc_content'].includes(r.id));
             s.adv.circular_track_slots_enabled=true;
             const r=s.adv.circular_track_slots.find(r=>r.id==='gc_content');r.width=value;r.radius=value;
             await window.Vue.nextTick();});
            }""", {'slots': seed['config']['adv']['circular_track_slots'], 'value': sample['input']})
            save_roundtrip(browser, page, sample['name'], True, {'width': sample['input'], 'radius': sample['input']}, sample['canonical'])
            page.close()
        # Rendered protection uses real Gallery artifacts and a schema-valid warning.
        page = fresh(browser)
        load(page, root / fixture)
        page.evaluate("""async ()=>{const s=(await import('./js/state.js')).state;
         s.annotationWarnings.value=[{code:'feature_selector_unmatched',setId:'probe',annotationId:'unmatched',recordId:'NC_001879.2',recordIndex:0,missingCount:1,message:'Skipped annotation: 1 feature selector(s) unmatched.',resultIndex:0,resultName:s.results.value[0].name}];
         await window.Vue.nextTick();}""")
        rendered = save_roundtrip(browser, page, 'rendered-metadata', False)
        assert rendered['runMetadata']['annotationWarnings'] and rendered['runMetadata']['trackSlotGeometry']['records']
        page.get_by_role('button', name='Linear', exact=True).click()
        page.wait_for_function("window.__GBDRAW_APP__.mode==='linear'")
        save_roundtrip(browser, page, 'inactive-biological', False)
        page.close()
        # Current writer admission rejects invalid/unfinished drafts, without null repair.
        page = fresh(browser)
        load(page, root / fixture)
        for name in ('invalid_text', 'incomplete_text', 'typed_invalid', 'typed_incomplete', 'typed_empty'):
            sample = next(s for s in samples if s['name'] == name)
            page.evaluate("""async value=>{const s=(await import('./js/state.js')).state;const r=s.adv.circular_track_slots.find(r=>r.id==='gc_content');r.width=value;r.radius=value;await window.Vue.nextTick();}""", sample['input'])
            before = page.evaluate(read)
            outcome = page.evaluate("""async ()=>{try {const r=await (await import('./js/services/config.js')).exportSession('invalid repair probe');return {accepted:r.status==='saved'};}catch(e){return {accepted:false,error:e.message};}}""")
            assert not outcome['accepted'] and page.evaluate(read) == before
            candidate = copy.deepcopy(rendered)
            row = next(r for r in candidate['config']['adv']['circular_track_slots'] if r['id'] == 'gc_content')
            row['width'] = row['radius'] = sample['input']
            path = args.artifacts / f'{name}.invalid.json'
            path.write_text(json.dumps(candidate))
            reader = fresh(browser)
            load(reader, root / fixture)
            prior = reader.evaluate(read)
            admission = load(reader, path, accepted=False)
            assert reader.evaluate(read) == prior
            reader.close()
            observations.append({'name': name, 'save': outcome, 'load': admission, 'rawDraft': sample['input'], 'saveStateUnchanged': True, 'failedLoadStateUnchanged': True})
        page.close()
        # Source-free invalid draft must also stay rejected at the changed boundary.
        page = fresh(browser)
        page.evaluate("""async ()=>{const s=(await import('./js/state.js')).state;s.adv.circular_track_slots_enabled=true;s.adv.circular_track_slots[0].width={value:'1e',unit:'px'};await window.Vue.nextTick();}""")
        before = page.evaluate(read)
        outcome = page.evaluate("""async ()=>{try {await (await import('./js/services/config.js')).exportSession('invalid settings probe');return {accepted:true};}catch(e){return {accepted:false,error:e.message};}}""")
        assert not outcome['accepted'] and page.evaluate(read) == before
        observations.append({'name': 'settings-invalid', 'save': outcome, 'saveStateUnchanged': True})
        page.close()
        browser.close()
finally:
    server.shutdown()
    server.server_close()
    thread.join()

assert not errors and not external, (errors, external)
result = {'sourceBaseSha': subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
          'sourcePatchSha256': sha(subprocess.check_output(['git', '-C', str(repo), 'diff', '--', 'gbdraw/web/js/services/config.js'])),
          'sourceSha256': {path: sha((root / path).read_bytes()) for path in ('gbdraw/web/js/services/config.js', 'gbdraw/web/js/services/session-authority.js', 'gbdraw/web/js/services/session-active-config-contract.js', 'gbdraw/web/js/app/circular-track-slots.js')},
          'snapshotMatchesCheckout': True, 'wheelPythonFilesMatched': len(python_paths), 'wheelSha256': sha(wheel.read_bytes()),
          'fixture': fixture, 'fixtureSha256': sha((root / fixture).read_bytes()), 'scalarFixtureSha256': sha(Path(__file__).with_name('scalar-fixtures.json').read_bytes()),
          'environment': {'python': platform.python_version(), 'playwright': version('playwright'), 'chromium': browser_version, 'nativePackage': str(Path(__import__('gbdraw').__file__).resolve()), 'url': url, 'sourceRoot': str(root)},
          'observations': observations, 'pageErrors': errors, 'externalRequests': external,
          'limits': ['No S01 codec, new controls, IME/unit/Auto lifetime, geometry change or native render replay.',
                     'Rendered warning is a synthetic valid metadata fixture; resolved geometry and SVG are existing public Gallery artifacts.',
                     'Scalar fixtures are assigned to existing reactive drafts; this does not test future DOM controls.',
                     'Native reader acceptance is separate from Web admission; invalid native config domain is unchanged.']}
args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
print('WROTE', args.output, flush=True)
