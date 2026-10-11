#!/usr/bin/env python3
"""Print a campaign's progress from <task dir>/progress.json.

Usage: progress.py <task dir> [--gh] [--write]

Each item has a weight and a status. Credit per status: todo 0, in_progress 0.3,
in_ci 0.8, merged and done 1; dropped items are left out. --gh reads the state
of each item's `pr` from GitHub (merged -> merged, closed -> dropped, open with
the item still todo -> in_ci); --write saves those updates back to the file.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess

CREDIT = {'todo': 0.0, 'in_progress': 0.3, 'in_ci': 0.8, 'merged': 1.0, 'done': 1.0}


def pr_state(number):
    out = subprocess.run(['gh', 'pr', 'view', str(number), '--json', 'state'], capture_output=True, text=True)
    return json.loads(out.stdout).get('state') if out.returncode == 0 else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('task_dir')
    parser.add_argument('--gh', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    path = pathlib.Path(args.task_dir) / 'progress.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    items = data.get('items', [])

    changed = False
    if args.gh:
        for item in items:
            if not item.get('pr') or item.get('status') in ('merged', 'done', 'dropped'):
                continue
            state = pr_state(item['pr'])
            new = {'MERGED': 'merged', 'CLOSED': 'dropped'}.get(state)
            if state == 'OPEN' and item.get('status') == 'todo':
                new = 'in_ci'
            if new and new != item.get('status'):
                item['status'] = new
                changed = True

    live = [item for item in items if item.get('status') != 'dropped']
    total = sum(float(item.get('weight', 1)) for item in live) or 1.0
    earned = sum(float(item.get('weight', 1)) * CREDIT.get(item.get('status', 'todo'), 0.0) for item in live)
    print(f"{data.get('campaign', path.parent.name)}: {100 * earned / total:.0f}% "
          f"({sum(1 for i in live if i.get('status') in ('merged', 'done'))} of {len(live)} items finished)")
    for item in items:
        pr = f" #{item['pr']}" if item.get('pr') else ''
        print(f"  [{item.get('status', 'todo'):>11}] w{item.get('weight', 1)} {item.get('id', '?')}{pr} {item.get('title', '')}")
    if changed and args.write:
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        print(f'updated {path}')


if __name__ == '__main__':
    main()
