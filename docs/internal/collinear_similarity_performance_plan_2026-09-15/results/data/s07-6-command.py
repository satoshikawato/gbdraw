"""Record one S07.6 command, full output, wall time, exit and host snapshots."""
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time

out = Path(__file__).resolve().parent
name, *command = sys.argv[1:]
logs = out / 's07-6-logs'
logs.mkdir(exist_ok=True)

def host():
    return subprocess.check_output(['ps', '-eo', 'pid,ppid,pcpu,etimes,comm', '--sort=-pcpu'], text=True)

start = datetime.datetime.now(datetime.timezone.utc).isoformat()
(logs / f'{name}-host-before.txt').write_text(host())
tick = time.monotonic()
with (logs / f'{name}.log').open('w') as log:
    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
entry = {'name': name, 'argv': command, 'cwd': str(Path.cwd()), 'startUtc': start,
         'endUtc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'wallSeconds': time.monotonic()-tick, 'exitCode': result.returncode,
         'affinity': sorted(os.sched_getaffinity(0)), 'PYTHONHASHSEED': os.environ.get('PYTHONHASHSEED'),
         'log': str((logs / f'{name}.log').relative_to(out))}
(logs / f'{name}-host-after.txt').write_text(host())
with (out / 's07-6-commands.jsonl').open('a') as stream:
    stream.write(json.dumps(entry)+'\n')
print(json.dumps(entry), flush=True)
raise SystemExit(result.returncode)
