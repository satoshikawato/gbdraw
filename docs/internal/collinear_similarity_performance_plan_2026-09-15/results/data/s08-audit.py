"""Hash final review categories and verify predecessor/shared-tree preservation."""
from collections import Counter
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

DATA = Path(__file__).resolve().parent
ROOT = DATA.parents[4]
SHARED = ROOT.parent.parent
PREVIOUS = ROOT.parent / 'collinear-s07-8-20260916'
OUTPUT = DATA / 's08-review.json'


def git(root, *args):
    return subprocess.check_output(['git', *args], cwd=root,
                                   env={**os.environ, 'GIT_OPTIONAL_LOCKS': '0'})


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def category(path):
    if '/results/data/' in path:
        return 'generated evidence and reproduction commands'
    if path.startswith('gbdraw/'):
        return 'production'
    if path.startswith('tests/'):
        return 'tests'
    if path.startswith('tools/'):
        return 'benchmark tooling'
    return 'documentation'


start = json.loads((DATA / 's08-start.json').read_text())
preservation = {}
for label, root, expected_head, status_key, diff_key in (
    ('shared', SHARED, '51d2ae362ac154a3e360727ad086ad1d27bbc989', 'sharedStatus', 'sharedTrackedDiffSha256'),
    ('predecessor', PREVIOUS, start['head'], 'previousStatus', 'previousTrackedDiffSha256'),
):
    head = git(root, 'rev-parse', 'HEAD').decode().strip()
    status = git(root, 'status', '--short').decode()
    patch = sha(git(root, 'diff', '--binary', 'HEAD'))
    assert head == expected_head, label
    assert status == start[status_key], label
    assert patch == start[diff_key], label
    preservation[label] = {
        'head': head, 'statusUnchanged': True, 'trackedDiffSha256': patch,
        'branch': git(root, 'branch', '--show-current').decode().strip(),
        'upstream': git(root, 'for-each-ref', '--format=%(upstream)', 'refs/heads/' +
                        git(root, 'branch', '--show-current').decode().strip()).decode().strip() or None,
    }

allowed_inherited_edits = {
    'docs/internal/collinear_similarity_performance_plan_2026-09-15/MASTER_PLAN.md',
    'tests/run_protein_local_browser_acceptance.py',
}
for entry in start['inheritedDirtyFiles']:
    assert sha((PREVIOUS / entry['path']).read_bytes()) == entry['sha256'], entry['path']
    if entry['path'] not in allowed_inherited_edits:
        assert sha((ROOT / entry['path']).read_bytes()) == entry['sha256'], entry['path']

branch = git(ROOT, 'branch', '--show-current').decode().strip()
assert branch == start['branch']
assert git(ROOT, 'rev-parse', 'HEAD').decode().strip() == start['head']
assert git(ROOT, 'rev-parse', 'origin/dev').decode().strip() == start['base']
assert not git(ROOT, 'for-each-ref', '--format=%(upstream)', 'refs/heads/' + branch).strip()
assert git(ROOT, 'merge-base', '--is-ancestor', start['base'], 'HEAD') == b''
assert not git(ROOT, 'diff', '--name-only', 'origin/dev', '--', 'tests/reference_outputs', 'examples/gbdraw_social_preview.png')
assert not git(ROOT, 'diff', '--check')

paths = set(git(ROOT, 'diff', '--name-only', 'origin/dev').decode().splitlines())
paths.update(git(ROOT, 'ls-files', '--others', '--exclude-standard').decode().splitlines())
paths.discard(str(OUTPUT.relative_to(ROOT)))
files = []
for path in sorted(paths):
    source = ROOT / path
    if source.is_file():
        files.append({'path': path, 'category': category(path),
                      'sha256': sha(source.read_bytes()), 'bytes': source.stat().st_size})
report = {
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'base': start['base'], 'head': start['head'], 'branch': branch, 'upstream': None,
    'preservation': preservation,
    'predecessorDirtyFilesUnchanged': len(start['inheritedDirtyFiles']),
    'inheritedFilesChangedByS08': sorted(allowed_inherited_edits),
    's08ProductionFiles': ['gbdraw/web/js/services/config.js'],
    'sourceIsHeadPlusReviewedWorkingTree': True,
    'referenceAndOwnerShowcaseUnchanged': True,
    'categoryCounts': dict(Counter(row['category'] for row in files)),
    'files': files,
    'selfExcludedFromHashManifest': str(OUTPUT.relative_to(ROOT)),
}
OUTPUT.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'files'}, indent=2))
