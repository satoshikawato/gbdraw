#!/usr/bin/env python3
"""SessionStart hook (compact, resume): print the bound campaign's STATE.md.

The gbdraw-campaign skill binds a session to a task folder by writing the
folder path to ~/.claude/campaign-sessions/<session id>. After a compaction
or a resume, this hook prints the first 150 lines of that folder's STATE.md,
so the current state comes from the file rather than from a summary.
Sessions without a binding get no output.
"""
from __future__ import annotations

import json
import pathlib
import sys

MAX_LINES = 150


def main():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return
    session = str(payload.get('session_id') or '')
    if not session or '/' in session:
        return
    binding = pathlib.Path.home() / '.claude' / 'campaign-sessions' / session
    try:
        state = pathlib.Path(binding.read_text(encoding='utf-8').strip()) / 'STATE.md'
        lines = state.read_text(encoding='utf-8').splitlines()
    except OSError:
        return
    print(f'Campaign state reloaded from {state} (first {min(len(lines), MAX_LINES)} of {len(lines)} lines):')
    print('\n'.join(lines[:MAX_LINES]))


if __name__ == '__main__':
    main()
