#!/usr/bin/env python3
"""PreToolUse hook (Bash) for gbdraw subagents: deny a Playwright run whose
output would stream into the agent's context.

A command passes when it goes through .claude/scripts/run_quiet.py, redirects
its output to a file, or pipes it into tail, head, or grep. Any other
`playwright test` or `npm run test:web...` command is denied with the reason.
"""
from __future__ import annotations

import json
import re
import sys

RUNNER = re.compile(r'\b(?:npx\s+playwright\s+test|playwright\s+test|npm\s+run\s+test:web)')
BOUNDED = re.compile(r'run_quiet\.py|(?<![&\d])>\s*[^&\s]|\d>\s*[^&\s]|\|\s*(?:tail|head|grep)\b')


def main():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return
    command = (payload.get('tool_input') or {}).get('command') or ''
    if not RUNNER.search(command) or BOUNDED.search(command):
        return
    print(json.dumps({
        'hookSpecificOutput': {
            'hookEventName': 'PreToolUse',
            'permissionDecision': 'deny',
            'permissionDecisionReason': (
                'Run Playwright through .claude/scripts/run_quiet.py '
                '(python3 .claude/scripts/run_quiet.py -- <command>) or redirect its output '
                'to a log file, so the full output stays out of the context.'
            ),
        }
    }))


if __name__ == '__main__':
    main()
