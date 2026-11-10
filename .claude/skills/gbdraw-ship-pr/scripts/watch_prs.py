#!/usr/bin/env python3
"""Print one line per pull-request event: merged, closed, merge conflict, failed
check, a check that ended cancelled once nothing is pending, or all checks green.

Usage: watch_prs.py --prs 940 941 [--prs-file FILE] [--repo OWNER/NAME] [--interval 300]
                    [--once] [--seen-file FILE]

Run it under Monitor or in the background; each printed line is one event.
`--prs-file` is re-read every cycle, so PR numbers can be added while it runs.
A head counts as green only when every workflow run of that head has finished: a
queued run has no check runs yet, so its checks alone would look green too early.
`--once` keeps no memory between calls; with `--seen-file FILE` the lines already
printed are kept in FILE, so a repeated `--once` prints only new events.
"""
from __future__ import annotations

import argparse
import collections
import json
import subprocess
import time

FAILED = ('failure', 'timed_out', 'action_required')


def gh_api(path):
    out = subprocess.run(['gh', 'api', path], capture_output=True, text=True)
    return json.loads(out.stdout) if out.returncode == 0 else None


def default_repo():
    out = subprocess.run(['gh', 'repo', 'view', '--json', 'nameWithOwner', '-q', '.nameWithOwner'],
                         capture_output=True, text=True)
    return out.stdout.strip() or 'satoshikawato/gbdraw'


def pr_numbers(args):
    numbers = list(args.prs)
    if args.prs_file:
        try:
            with open(args.prs_file, encoding='utf-8') as handle:
                numbers += handle.read().split()
        except OSError:
            pass
    return list(dict.fromkeys(numbers))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prs', nargs='*', default=[])
    parser.add_argument('--prs-file')
    parser.add_argument('--repo')
    parser.add_argument('--interval', type=int, default=300, help='seconds between polls (default 300)')
    parser.add_argument('--once', action='store_true', help='poll once and exit')
    parser.add_argument('--seen-file', help='keep the printed event lines here, so a later run prints only new ones')
    args = parser.parse_args()
    repo = args.repo or default_repo()
    seen = set()
    if args.seen_file:
        try:
            with open(args.seen_file, encoding='utf-8') as handle:
                seen.update(handle.read().splitlines())
        except OSError:
            pass

    def emit(key, line):
        if key in seen or line in seen:
            return
        seen.add(key)
        seen.add(line)
        print(line, flush=True)
        if args.seen_file:
            with open(args.seen_file, 'a', encoding='utf-8') as handle:
                handle.write(line + '\n')

    while True:
        for n in pr_numbers(args):
            pr = gh_api(f'repos/{repo}/pulls/{n}')
            if not pr:
                continue
            sha = pr['head']['sha']
            if pr.get('merged'):
                emit((n, 'm'), f'PR #{n} MERGED')
                continue
            if pr['state'] == 'closed':
                emit((n, 'c'), f'PR #{n} CLOSED unmerged')
                continue
            if pr.get('mergeable_state') == 'dirty':
                emit((n, sha, 'd'), f'PR #{n} has merge conflicts (dirty)')
            runs = gh_api(f'repos/{repo}/commits/{sha}/check-runs?per_page=100')
            if not runs:
                continue
            by_name = collections.defaultdict(list)
            for run in runs['check_runs']:
                if '${{' in run['name']:  # unexpanded matrix placeholder of a cancelled run
                    continue
                by_name[run['name']].append(run)
            latest = {name: max(rs, key=lambda r: r['id']) for name, rs in by_name.items()}
            workflow_runs = gh_api(f'repos/{repo}/actions/runs?head_sha={sha}&per_page=100')
            runs_open = workflow_runs is None or any(r['status'] != 'completed' for r in workflow_runs['workflow_runs'])
            pending = runs_open or any(r['status'] != 'completed' for r in latest.values())
            fails = sorted(name for name, r in latest.items() if r['conclusion'] in FAILED)
            if fails:
                emit((n, sha, 'f', tuple(fails)), f'PR #{n} ({sha[:8]}) FAILED: {", ".join(fails)}')
            if not pending:
                cancelled = sorted(name for name, r in latest.items() if r['conclusion'] == 'cancelled')
                if cancelled:
                    emit((n, sha, 'x'), f'PR #{n} ({sha[:8]}) CANCELLED without success: {", ".join(cancelled)}')
                elif not fails:
                    emit((n, sha, 'g'), f'PR #{n} ({sha[:8]}) all checks completed green')
        if args.once:
            return
        time.sleep(args.interval)


if __name__ == '__main__':
    main()
