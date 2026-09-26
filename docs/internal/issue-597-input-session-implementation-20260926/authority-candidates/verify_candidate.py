"""Verify a candidate in a clean dev worktree or on a temporary policy copy."""

from copy import deepcopy
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def require(condition, message):
    if not condition:
        raise SystemExit(message)


inert = '--inert' in sys.argv[1:]
mode, worktree_arg = [arg for arg in sys.argv[1:] if arg != '--inert']
require(not inert or mode == 'privileged', '--inert supports the permission candidate only')
require(mode in ('product', 'privileged'), 'Use: verify_candidate.py product|privileged <dev-worktree>')
worktree = Path(worktree_arg).resolve()
plan = Path(__file__).resolve().parent.parent


def git(*args):
    return subprocess.check_output(['git', '-C', str(worktree), *args], text=True)


target = ('docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'
          if mode == 'product' else 'tools/web-change-policy.json')
if not inert:
    require(git('diff', '--name-only', 'HEAD').splitlines() == [target], 'Unexpected changed paths')
    require(not git('ls-files', '--others', '--exclude-standard').strip(), 'Unexpected untracked files')
base = git('show', f'HEAD:{target}')
if inert:
    patch = plan / 'authority-candidates/02-import-worker-permission.patch'
    subprocess.run(['git', 'apply', '--check', str(patch)], cwd=worktree, check=True)
    with tempfile.TemporaryDirectory(prefix='issue597-permission-copy-') as directory:
        copy = Path(directory) / target
        copy.parent.mkdir(parents=True)
        copy.write_text(base)
        subprocess.run(['git', 'apply', str(patch)], cwd=directory, check=True)
        candidate = copy.read_text()
else:
    candidate = (worktree / target).read_text()

if mode == 'product':
    fields = {
        'Concern': 'concern', 'Scenario revision': 'scenarioRevision', 'Choice': 'choice',
        'Rationale': 'rationale', 'Must preserve': 'mustPreserve', 'May retire': 'mayRetire',
        'Accepted residual risk': 'acceptedResidualRisk', 'Owner': 'owner', 'Decision date': 'decisionDate',
    }
    record_pattern = r'^### (PD-OI-\d+):.*?(?=^### PD-OI-|^## Acceptance contract catalog|\Z)'
    def records(text):
        return {match.group(1): match.group(0) for match in re.finditer(record_pattern, text, re.S | re.M)}

    before, after = records(base), records(candidate)
    require(set(after) - set(before) == {'PD-OI-044', 'PD-OI-045'}, 'Wrong new decision IDs')
    require(all(after.get(key) == value for key, value in before.items()), 'Existing outcome changed')
    def clauses(text):
        return text.split('## Interpretation and lifecycle', 1)[1].split('## Product Decision records', 1)[0]
    require(clauses(base) == clauses(candidate), 'Existing lifecycle/cross-surface clauses changed')
    revision = int(re.search(r'- Contract revision: `(\d+)`', base).group(1))
    require(f'- Contract revision: `{revision + 1}`' in candidate, 'Wrong contract revision')
    require('diagram-generation.circular-source-collection' not in candidate, 'BUG-01 outcome included')
    map_data = json.loads(git('show', 'HEAD:tools/web-product-impact-map.json'))
    decision_authority = json.loads(git('show', 'HEAD:tools/web-product-decisions.json'))
    for record_id, filename in [('PD-OI-044', '02_RECORD_DISCOVERY.md'), ('PD-OI-045', '03_SESSION_OPERATIONS.md')]:
        source = (plan / 'decisions' / filename).read_text()
        human = re.search(r'```text\nPRODUCT_DECISION\n(.*?)\n```', source, re.S).group(1)
        receipt = {}
        for line in human.splitlines():
            label, value = line.split(': ', 1)
            receipt[fields[label]] = int(value) if label == 'Scenario revision' else value
        serialized = json.loads(re.search(r'```json\n(.*?)\n```', source, re.S).group(1))
        published = json.loads(re.search(r'```json\n(.*?)\n```', after[record_id], re.S).group(1))
        require(len(receipt) == 9 and receipt == serialized == published, f'{record_id}: field mismatch')
        for label, value in [('Concern key', receipt['concern']), ('Scenario revision', receipt['scenarioRevision']), ('Selected outcome', receipt['choice'])]:
            require(f'- {label}: `{value}`' in after[record_id], f'{record_id}: metadata mismatch')
        require(all(value for value in receipt.values()), f'{record_id}: empty receipt field')
        require(receipt['owner'] in decision_authority['maintainerLogins'], 'Owner not in base allowlist')
        require(not any(c['key'] == receipt['concern'] for c in map_data['concerns']), 'Concern now mapped; revisit target')
        require(not any(d['concernKey'] == receipt['concern'] for d in decision_authority['decisions']), 'Durable authority conflict')
        digest = re.search(r'承認対象資料の SHA-256: `([a-f0-9]{64})`', source).group(1)
        require(digest in after[record_id], 'Approval provenance digest missing')
        print(f'{record_id}: all 9 human/serialized/candidate fields equal; owner and provenance valid')
    marker = '## Acceptance contract catalog'
    require(base.split(marker, 1)[1] == candidate.split(marker, 1)[1], 'Existing acceptance/risk text changed')
else:
    before, after = json.loads(base), json.loads(candidate)
    expected = deepcopy(before)
    expected['allowedPrivilegedOwners']['Diagram Worker'].append('services/session-import-client.js')
    expected['allowedPrivilegedOwners']['Diagram Worker'].sort()
    expected['allowedPrivilegedImporters']['services/session-file.js'].append('workers/session-import-worker.js')
    expected['allowedPrivilegedImporters']['services/session-file.js'].sort()
    require(after == expected, 'Permission exceeds one constructor owner and one codec import edge')
    # These source strings are detector probes; they never construct a Worker.
    probe = r'''
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { detectPrivilegedWebCapabilities } from './tools/web-architecture-detectors.mjs';
const sources = {
  'services/config.js': "import { readSessionFile } from './session-import-client.js';",
  'services/session-import-client.js': "export const readSessionFile = () => new Worker(new URL('../workers/session-import-worker.js', import.meta.url), { type: 'module' });",
  'workers/session-import-worker.js': "import { readSessionText } from '../services/session-file.js'; self.onmessage = async ({ data }) => self.postMessage(JSON.parse(await readSessionText(data.file)));"
};
const detected = detectPrivilegedWebCapabilities(sources);
assert.deepEqual(detected.operatorMatchesByCapability['Diagram Worker'], [{ path: 'services/session-import-client.js', count: 1 }]);
for (const [name, matches] of Object.entries(detected.operatorMatchesByCapability)) {
  if (name !== 'Diagram Worker') assert.deepEqual(matches, []);
}
assert.deepEqual(detected.importersByTarget['services/session-file.js'], ['workers/session-import-worker.js']);
assert.equal(Object.values(detected.importersByTarget).flat().length, 1);
const base = JSON.parse(execFileSync('git', ['show', 'HEAD:tools/web-change-policy.json'], { encoding: 'utf8' }));
const candidate = JSON.parse(readFileSync('tools/web-change-policy.json', 'utf8'));
assert.equal(base.allowedPrivilegedOwners['Diagram Worker'].includes('services/session-import-client.js'), false);
assert.equal(candidate.allowedPrivilegedOwners['Diagram Worker'].includes('services/session-import-client.js'), true);
assert.equal(candidate.allowedPrivilegedOwners['Diagram Worker'].includes('services/other-worker-client.js'), false);
assert.equal(base.allowedPrivilegedImporters['services/session-file.js'].includes('workers/session-import-worker.js'), false);
assert.equal(candidate.allowedPrivilegedImporters['services/session-file.js'].includes('workers/session-import-worker.js'), true);
console.log('Detector probe: one constructor owner; one codec import edge; unrelated owner denied');
'''
    if inert:
        probe = probe.replace("JSON.parse(readFileSync('tools/web-change-policy.json', 'utf8'))", json.dumps(after))
    subprocess.run(['node', '--input-type=module'], input=probe, text=True, cwd=worktree, check=True)
print(f'{mode}: candidate shape and single target path PASS')
