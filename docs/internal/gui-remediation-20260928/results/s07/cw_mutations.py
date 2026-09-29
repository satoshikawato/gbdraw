"""Apply one temporary CW mutation, run the CW spec, restore, and verify SHA-256."""
import hashlib
import os
import pathlib
import subprocess
import sys

WT = pathlib.Path(__file__).resolve().parents[5]
MUTATIONS = {
    # M2 (CW-01): a status computed calls the canonical request builder again.
    'M2': ('gbdraw/web/js/app/app-setup.js', [
        ("import { createRulePreparation } from './rule-matching.js';\n",
         "import { createRulePreparation } from './rule-matching.js';\n"
         "import { buildCanonicalRenderRequest as cwMutationBuild } from '../services/session-request.js';\n"),
        ("  const linearComparisonUi = computed(() => projectLinearComparisonUi({\n",
         "  const linearComparisonUi = computed(() => { try { cwMutationBuild({ state, filesData: {} }); } catch (_error) {} return projectLinearComparisonUi({\n"),
    ]),
    # M3 (CW-02): the request builder ignores the label table built for this Generate.
    'M3': ('gbdraw/web/js/services/session-request.js', [
        ("    || (typeof state.generatedLabelOverrideTsv === 'string' ? state.generatedLabelOverrideTsv",
         "    || (false ? state.generatedLabelOverrideTsv"),
    ]),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(name):
    relative, edits = MUTATIONS[name]
    path = WT / relative
    original = path.read_bytes()
    before = sha(path)
    text = original.decode()
    for old, new in edits:
        assert text.count(old) == 1, (name, old)
        text = text.replace(old, new, 1)
    if name == 'M2':
        # Close the block body opened by the mutation.
        anchor = text.index('cwMutationBuild({ state, filesData: {} })')
        end = text.index('}));\n', anchor)
        text = text[:end] + '}); });\n' + text[end + len('}));\n'):]
    path.write_text(text)
    try:
        env = dict(os.environ, NODE_PATH=os.environ.get('NODE_PATH', ''))
        run = subprocess.run(['npx', 'playwright', 'test', 'tests/web/computation-ownership.playwright.spec.js',
                              '--workers=1', '--reporter=line'], cwd=WT, env=env, capture_output=True, text=True,
                             timeout=1800)
        lines = [line for line in (run.stdout + run.stderr).splitlines()
                 if 'passed' in line or 'failed' in line or 'Expected' in line or 'Received' in line or 'Error:' in line]
        print(name, 'exit', run.returncode)
        print('\n'.join(lines[:14]))
    finally:
        path.write_bytes(original)
        print(name, 'restored', sha(path) == before, before[:16])


if __name__ == '__main__':
    main(sys.argv[1])
