"""Record exact verification commands; elapsed time is not performance evidence."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[5]
name, *argv = sys.argv[1:]
row = {'name': name, 'argv': argv, 'cwd': str(ROOT),
       'start': datetime.now(timezone.utc).isoformat()}
with (OUT / f'{name}.log').open('w') as log:
    result = subprocess.run(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
row.update(end=datetime.now(timezone.utc).isoformat(), exit=result.returncode)
with (OUT / 'commands.jsonl').open('a') as log:
    log.write(json.dumps(row) + '\n')
print(json.dumps(row), flush=True)
raise SystemExit(result.returncode)
