#!/usr/bin/env python3
"""Bind the current Claude Code session to a campaign task folder.

Usage: bind_session.py <task dir>            bind $CLAUDE_CODE_SESSION_ID
       bind_session.py --clear <task dir>    remove every binding to the folder

The binding is ~/.claude/campaign-sessions/<session id>, holding the folder's
absolute path. The SessionStart hook .claude/hooks/campaign-state.py reads it
and prints the folder's STATE.md after a compaction or a resume.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys

BINDINGS = pathlib.Path.home() / '.claude' / 'campaign-sessions'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('task_dir')
    parser.add_argument('--clear', action='store_true')
    args = parser.parse_args()
    folder = pathlib.Path(args.task_dir).expanduser().resolve()

    if args.clear:
        removed = 0
        for binding in BINDINGS.glob('*') if BINDINGS.is_dir() else []:
            if binding.read_text(encoding='utf-8').strip() == str(folder):
                binding.unlink()
                removed += 1
        print(f'removed {removed} binding(s) to {folder}')
        return

    session = os.environ.get('CLAUDE_CODE_SESSION_ID', '')
    if not session or '/' in session:
        sys.exit('CLAUDE_CODE_SESSION_ID is not set; run this from a Claude Code Bash tool call.')
    if not (folder / 'STATE.md').is_file():
        print(f'warning: {folder}/STATE.md does not exist yet', file=sys.stderr)
    BINDINGS.mkdir(parents=True, exist_ok=True)
    (BINDINGS / session).write_text(f'{folder}\n', encoding='utf-8')
    print(f'session {session} -> {folder}')


if __name__ == '__main__':
    main()
