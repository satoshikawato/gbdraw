import datetime, hashlib, json, os, pathlib, subprocess, sys, time
root = pathlib.Path('/tmp/issue597-S05-trusted-guard-integration-20260927-evidence')
root.mkdir(exist_ok=True)
name = sys.argv[1]
cmd = sys.argv[2:]
cwd = '/tmp/gbdraw-issue597-S05.ujOWyl'
log = root / (name + '.log')
if log.exists():
    raise SystemExit('Refusing to overwrite invocation: ' + str(log))
start = datetime.datetime.now(datetime.timezone.utc).isoformat()
t0 = time.monotonic()
with log.open('wb') as stream:
    result = subprocess.run(cmd, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT)
record = {'name': name, 'command': cmd, 'cwd': cwd, 'startedAt': start,
          'completedAt': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'durationSeconds': time.monotonic()-t0, 'exit': result.returncode,
          'log': str(log), 'sha256': hashlib.sha256(log.read_bytes()).hexdigest()}
(root / (name + '.json')).write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
print(log.read_text(errors='replace')[-5500:])
sys.exit(result.returncode)
