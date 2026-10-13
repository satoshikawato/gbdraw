#!/bin/bash
# Wait once until a sub-agent of another live session hands back (its transcript ends with a SubagentHandback call).
# Usage: wait_handback.sh LABEL TRANSCRIPT.jsonl PARENT_PID REPORT_FILE [MAX_MINUTES]
label=$1; f=$2; pid=$3; report=$4; max=${5:-150}; start=$(date +%s)
while :; do
  last=$(tail -n 6 "$f" 2>/dev/null | /home/kawato/micromamba/bin/python3.13 -c '
import json, sys
name = ""
for line in sys.stdin:
    try: c = (json.loads(line).get("message") or {}).get("content")
    except Exception: continue
    if isinstance(c, list):
        for x in c:
            if x.get("type") == "tool_use": name = x.get("name") or ""
print(name)')
  if [ "$last" = SubagentHandback ]; then echo "$label: agent handed back at $(date +%T)"; break; fi
  if ! kill -0 "$pid" 2>/dev/null; then echo "$label: parent process $pid is gone, no hand-back (last transcript write $(stat -c %y "$f" | cut -c12-19))"; break; fi
  if [ $(( $(date +%s) - start )) -gt $(( max * 60 )) ]; then echo "$label: still running after $max minutes (last transcript write $(stat -c %y "$f" | cut -c12-19))"; break; fi
  sleep 60
done
ls -la --time-style=+%H:%M:%S "$report" 2>&1
