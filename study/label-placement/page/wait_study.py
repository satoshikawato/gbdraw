"""Wait once for the seven study sub-agents (parent session 5544d847) to finish or stall.

The lead agent a2b10cfa stopped PARTIAL at 00:19 on 2026-10-12; these seven write the designs, the sweep and the
survey. An agent counts as finished when its last transcript entry is a text-only assistant message and the
transcript has been idle 3 min, as stalled when it has been idle 30 min. Exits when none is still running
(or after 4 h) and prints one line per agent with its report file.
"""
import json, os, time

T = '/home/kawato/gbdraw-baselines/label-placement-20261011'
D = '/home/kawato/.claude/projects/-home-kawato-gbdraw-work/5544d847-7c4d-4337-99d7-e3d43f7fa157/subagents/'
AGENTS = {
    'A': ('agent-aa7761bfcb1e4acfc', 'logs/design-A.md'),
    'B': ('agent-a4d987edf076862c8', 'logs/design-B.md'),
    'C': ('agent-a24c401c6a3fdea98', 'logs/design-C.md'),
    'D': ('agent-ab05cfd1555ec4b54', 'logs/design-D.md'),
    'E': ('agent-a9c2b31bfee210c01', 'logs/design-E.md'),
    'sweep': ('agent-a99417282348b9f4d', 'logs/part2-breaking-points.md'),
    'survey': ('agent-a063181e53e2374c7', 'logs/survey-tools.md'),
}


def state(name):
    tr = D + name + '.jsonl'
    idle = time.time() - os.path.getmtime(tr)
    kinds = []
    try:
        with open(tr, 'rb') as f:
            f.seek(max(0, os.path.getsize(tr) - 1500000))
            last = json.loads(f.read().splitlines()[-1])
        content = (last.get('message') or {}).get('content')
        if last.get('type') == 'assistant' and isinstance(content, list):
            kinds = [c.get('type') for c in content]
    except Exception:
        pass
    if 'text' in kinds and 'tool_use' not in kinds and idle > 180:
        return 'finished', idle
    if idle > 1800:
        return 'stalled', idle
    return 'running', idle


start = time.time()
while True:
    states = {k: state(v[0]) for k, v in AGENTS.items()}
    timeout = time.time() - start > 4 * 3600
    if timeout or all(s != 'running' for s, _ in states.values()):
        print('TIMEOUT' if timeout else 'ALL SETTLED', f'waited={(time.time() - start) / 60:.0f}min')
        for k, (s, idle) in states.items():
            report = AGENTS[k][1]
            print(f'{k}: {s} idle={idle / 60:.0f}min report={"yes" if os.path.exists(T + "/" + report) else "MISSING"} ({report})')
        break
    time.sleep(60)
