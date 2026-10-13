#!/bin/bash
# Wait once until FILE exists (written by another session), the other session's process is gone, or MAX minutes pass.
# Usage: wait_file.sh FILE PARENT_PID [MAX_MINUTES]
f=$1; pid=$2; max=${3:-150}; start=$(date +%s)
until [ -s "$f" ]; do
  if ! kill -0 "$pid" 2>/dev/null; then echo "process $pid is gone; $f was not written"; exit 0; fi
  if [ $(( $(date +%s) - start )) -gt $(( max * 60 )) ]; then echo "$f not written after $max minutes"; exit 0; fi
  sleep 30
done
sleep 5; ls -la --time-style=+%H:%M:%S "$f"
