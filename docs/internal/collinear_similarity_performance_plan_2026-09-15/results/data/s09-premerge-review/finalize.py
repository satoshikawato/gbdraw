"""Finish documentation only after the S09 source/preservation audit passes."""
import difflib
import hashlib
import json
from pathlib import Path
import subprocess

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[5]
PLAN = OUT.parent.parent.parent
audit = json.loads((OUT / 'audit.json').read_text())
inheritance = json.loads((OUT / 'inheritance.json').read_text())
for name, expected in (audit['sourceSha256'] | audit['testsSha256']).items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
result = PLAN / 'results/S09_PREMERGE_REVIEW.md'
text = result.read_text().replace(
    'Final evidence audit is pending.',
    'Review, focused verification and final preservation audit are complete.')
text = text.replace('Final file-preservation audit is pending at draft.',
    'Final audit verifies all 208 wheel Python files, 466 installed package files,\n'
    '  inherited evidence, and unchanged HEAD/status/tracked diffs in 11 protected\n'
    '  worktrees. Current source/test hashes and separate diff categories are recorded.')
result.write_text(text)
master = PLAN / 'MASTER_PLAN.md'
master.write_text(master.read_text().replace('最終監査は実行中。',
    '最終source・wheel/package・継承証拠・既存11 worktreeの保持監査も完了した。'))
parts = []
documentation = audit['s09Diff']['documentation']
for name in documentation:
    original = Path(inheritance['source']) / name
    before = original.read_text() if original.exists() else ''
    parts.extend(difflib.unified_diff(before.splitlines(True), (ROOT / name).read_text().splitlines(True),
                 fromfile='inherited-final/' + name, tofile='s09/' + name))
(OUT / 'new-documentation.patch').write_text(''.join(parts))
subprocess.run(['git', 'diff', '--check'], cwd=ROOT, check=True)
completion = {'sourceAggregateSha256': audit['sourceAggregateSha256'],
    'auditSha256': hashlib.sha256((OUT / 'audit.json').read_bytes()).hexdigest(),
    'finalDocumentationSha256': {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest() for n in documentation},
    'sourceAndTestHashesStillMatchAudit': True,
    'humanReviewAndCommitFixing': 'Still required; no approval or publication performed.'}
(OUT / 'completion.json').write_text(json.dumps(completion, indent=2) + '\n')
print(json.dumps(completion, indent=2))
