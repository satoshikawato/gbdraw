"""Run retained SESSION 05A13-B programs with the permanent visual comparator.

The immutable archive is an input, never copied or modified. The original
statement runner keeps failed checkpoints and continues independent operations.
Nonzero exit status prevents its diagnostic continuation from becoming a PASS.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

# Retained programs are immutable inputs, including their archive directory.
sys.dont_write_bytecode = True

from playwright.sync_api import sync_playwright

HASHES = {
    'audit.py': 'd24c4627013cde32a57aa157c015788e1b4caabe5a7ac8265041f883bcea85f1',
    'harness.py': '96f49c44fe172d7eeec374bc2e11d19b174daee9199961844ec2aac2ca2efb51',
    'journeys.py': '0ec1adbac279ff7579985e416a269861e7006c6c2aef0c357bb635d84954eb07',
    'audit-state.js': 'e4525f9f4f56fdec788b8c3bb1141e1fa36da45488f22d00b8bc645b231b4cc8',
    'disposition-adapters.py': '30a225fba3ce9bc23eee10dcaf4c1c39ecb8150021bbd4cc3ba93679b78c98a2',
    'supplemental-controls.py': 'f40b4355817ea4248bc630587d1d4537b3cc86b49d18e40133653eee6a83c827',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--journeys', nargs='+', required=True)
    args = parser.parse_args()
    for name, expected in HASHES.items():
        assert hashlib.sha256((args.archive / name).read_bytes()).hexdigest() == expected, name
    for seed in json.loads((args.archive / 'seed-manifest.json').read_text()):
        assert hashlib.sha256(Path(seed['path']).read_bytes()).hexdigest() == seed['sha256'], seed['seed']
    args.output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.archive))
    import harness
    import journeys

    # Keep the original interpreter and all journey statements. Replace only its
    # comparator calls; both normal tests and this runner import the same helper.
    source = (args.archive / 'audit.py').read_text()
    for left, right in [('sm', 'mm'), ('bs', 'bm'), ('em', 'bs'), ('em', 'bm')]:
        # Binding enrichment is directional: selected -> mounted/export.
        expected, actual = (right, left) if left == 'em' else (left, right)
        source = source.replace(f'if {left}!={right}:', f'if not self.visual_equal({expected},{actual},d):')
    # Session 42 adds a source-free document, not a nullable canonical owner.
    # Retain the original writer checks for every source-bearing saved Session.
    old_oracle = "ok=schemas['session']==41 and schemas['request']==7 and schemas['bindings']==2 and (schemas['catalog']==3 or (cat is None and not d.get('results')))"
    assert old_oracle in source
    source = source.replace(old_oracle, 'ok=self.saved_session_valid(d,schemas)')
    spec = importlib.util.spec_from_file_location('retained_audit', args.archive / 'audit.py')
    audit = importlib.util.module_from_spec(spec)
    exec(compile(source, str(args.archive / 'audit.py'), 'exec'), audit.__dict__)
    state_source = audit.STATE
    start = state_source.index(' const props=')
    end = state_source.index(' if(contentOverride!==null)')
    audit.STATE = state_source[:start] + " const {svgVisualSemantics:semantics}=await import('/tests/web/helpers/svg-visual-semantics.mjs');\n" + state_source[end:]

    def visual_equal(self, left, right, state):
        catalog = state.get('catalog') or {}
        item = next((item for item in catalog.get('items', []) if item['resultIndex'] == state['selectedIndex']), {})
        return self.page.evaluate("""async ({left,right,ids})=>{
          const {compareVisualSemantics}=await import('/tests/web/helpers/svg-visual-semantics.mjs');
          return compareVisualSemantics(left,right,{catalogIds:ids}).length===0;
        }""", {'left': left, 'right': right, 'ids': [feature['svgId'] for feature in item.get('features', [])]})
    audit.AuditJourney.visual_equal = visual_equal

    def saved_session_valid(self, document, schemas):
        if schemas['session'] != 42 or schemas['bindings'] != 2:
            return False
        if document.get('renderRequest', 'missing') is not None:
            return schemas['request'] == 7 and (schemas['catalog'] == 3
                or (schemas['catalog'] is None and not document.get('results')))
        # Explicit settings, input/Result absence and honest committed-owner
        # absence supplement the shared document validator. Empty SVG equality
        # is never used as the settings-only oracle.
        return self.page.evaluate('''async document => {
          const {adoptCurrentSessionDocument,hasBiologicalSessionInputs}=await import('./js/services/session-authority.js');
          const {buildConfigData,getCommittedCanonicalSession}=await import('./js/services/config.js');
          const {state:s}=await import('./js/state.js');
          try {
            const admitted=adoptCurrentSessionDocument(document,42);
            return admitted.canonical===null && getCommittedCanonicalSession()===null
              && document.results.length===0 && s.results.value.length===0
              && document.editorState.featureCatalog===null && s.featureCatalog.value===null
              && !hasBiologicalSessionInputs({...s.files,linearSeqs:s.linearSeqs})
              && document.ui.mode===s.mode.value
              && JSON.stringify(document.config)===JSON.stringify(buildConfigData());
          } catch {return false;}
        }''', document)

    audit.AuditJourney.saved_session_valid = saved_session_valid
    original_edit = audit.AuditJourney.edit
    original_observe = audit.AuditJourney.observe

    def edit(self, kind, *params, **options):
        target = original_edit(self, kind, *params, **options)
        if self.jid == 'J13' and kind == 'placement':
            self.placement_target = target['svg_id']
        return target

    def observe(self, reason):
        state = original_observe(self, reason)
        if self.jid != 'J13' or not state or not hasattr(self, 'placement_target'):
            return state
        if reason.startswith(('before Undo ', 'before Redo ', 'after Undo ', 'after Redo ')):
            current = json.loads(gzip.decompress((self.folder / state['mounted_semantics']['path']).read_bytes()))
            if reason.startswith('before '):
                self.history_non_targets = current
            else:
                differences = self.page.evaluate("""async ({before,after,target})=>{
                  const {nonTargetFeatures,compareVisualSemantics}=await import('/tests/web/helpers/svg-visual-semantics.mjs');
                  return compareVisualSemantics(nonTargetFeatures(before,[target]),nonTargetFeatures(after,[target]));
                }""", {'before': self.history_non_targets, 'after': current, 'target': self.placement_target})
                self.record('non_target_preservation', 'OBSERVATION' if differences else 'PASS',
                            reason=reason, differences=differences[:5])
        return state

    audit.AuditJourney.edit = edit
    audit.AuditJourney.observe = observe
    harness.ROOT = journeys.ROOT = args.root.resolve()
    harness.OUT = audit.OUT = args.output.resolve()
    harness.LINEAR_FILES = journeys.LINEAR_FILES = [harness.ROOT / 'examples' / name for name in ['MellatMJNV.gb', 'MeenMJNV.gb', 'LvMJNV.gb']]
    # The retained harness uses port 4194 for its local-only request policy.
    harness.INIT = harness.INIT.replace('window.__SWEEP__=', "window.__DEFINITION_COMPLETED__=0;const log=console.log;console.log=(...a)=>{if(a[0]==='Definition text updated')window.__DEFINITION_COMPLETED__++;log(...a)};window.__SWEEP__=")
    original_load = audit.AuditJourney.load

    def load(self, seed):
        result = original_load(self, seed)
        if self.jid == 'C01' and '-03-save' in str(seed):
            self.page.wait_for_function('window.__DEFINITION_COMPLETED__ > 0', timeout=180000)
            self.observe('Web Load definition callback completed before Generate')
        return result
    audit.AuditJourney.load = load
    manifest = {'archive': str(args.archive.resolve()), 'sourceHashes': HASHES,
                'comparatorSha256': hashlib.sha256((args.root / 'tests/web/helpers/svg-visual-semantics.mjs').read_bytes()).hexdigest(), 'journeys': args.journeys}
    (args.output / 'runner-inputs.json').write_text(json.dumps(manifest, indent=2))
    failed = []
    with sync_playwright() as pw:
        for jid in args.journeys:
            if jid not in {f'J{i:02}' for i in range(1, 49)} | {'C01'}:
                raise ValueError(jid)
            browser = pw.chromium.launch(headless=True)
            journey = audit.AuditJourney(browser, jid, {'width': 390, 'height': 844} if jid in ['J43', 'J44'] else None)
            key = jid.lower()
            fn = audit.supplements.get(key) or audit.adapters.get(key) or getattr(journeys, key, None) or getattr(harness, key)
            # Original supplement OUT still references the immutable CLI fixtures.
            program = audit.Program(journey, fn)
            try:
                program.run()
                journey.observe('final')
                journey.snap('final', True)
                status = ('FAIL' if journey.observations else 'BLOCKED' if journey.blocked
                          else 'PASS_WITH_DISPOSITION' if journey.checkpoints else 'PASS')
            except Exception:
                import traceback
                journey.record('runner', 'HARNESS_ERROR', error=traceback.format_exc())
                status = 'HARNESS_ERROR'
            journey.finish(status)
            if status not in ['PASS', 'PASS_WITH_DISPOSITION']:
                failed.append(jid)
            browser.close()
    print(json.dumps({'failed': failed, 'completed': args.journeys}))
    return bool(failed)


if __name__ == '__main__':
    sys.exit(main())
