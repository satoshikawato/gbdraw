#!/usr/bin/env python3
"""Run a command with its full output in a log file and print a short summary.

Usage:
  run_quiet.py [--log FILE] [--tail N] -- <command> [args...]
  run_quiet.py [--log FILE] [--tail N] --shell "<shell command line>"

The summary shows the exit status, the duration, the log path, the lines that
look like failures or result counts (at most 40), and the last N lines
(default 15). The exit status of the command is preserved. Without --log, the
log goes to $GBDRAW_LOG_DIR, else $TMPDIR/run_quiet, else ./run_quiet-logs.
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import pathlib
import re
import subprocess
import sys
import time

SIGNAL = re.compile(
    r'(\bFAILED\b|\bERROR\b|\bError:|✘|\bnot ok\b|^# (?:fail|pass|tests)\b|'
    r'\b\d+ (?:passed|failed|flaky|skipped|errors?)\b|Traceback|AssertionError|Timeout|'
    r'^\s*\d+\) \[|^=+ .* in [\d.]+s)',
)


def default_log(slug):
    base = os.environ.get('GBDRAW_LOG_DIR') or (
        os.path.join(os.environ['TMPDIR'], 'run_quiet') if os.environ.get('TMPDIR') else 'run_quiet-logs')
    stamp = dt.datetime.now().strftime('%Y%m%d-%H%M%S')
    return pathlib.Path(base) / f'{stamp}-{slug}.log'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--log')
    parser.add_argument('--tail', type=int, default=15)
    parser.add_argument('--shell', help='run this command line through bash')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not args.shell and not command:
        parser.error('give a command after -- or use --shell')
    label = args.shell or ' '.join(command)
    slug = re.sub(r'[^A-Za-z0-9]+', '-', label)[:60].strip('-') or 'command'
    log = pathlib.Path(args.log) if args.log else default_log(slug)
    log.parent.mkdir(parents=True, exist_ok=True)

    start = time.monotonic()
    with log.open('w', encoding='utf-8', errors='replace') as handle:
        if args.shell:
            proc = subprocess.run(['bash', '-c', args.shell], stdout=handle, stderr=subprocess.STDOUT)
        else:
            proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    seconds = time.monotonic() - start

    lines = log.read_text(encoding='utf-8', errors='replace').splitlines()
    signal = [line for line in lines if SIGNAL.search(line)]
    print(f'$ {label}')
    print(f'exit {proc.returncode} after {seconds:.0f}s; {len(lines)} lines in {log}')
    if signal:
        print(f'--- matching lines ({min(len(signal), 40)} of {len(signal)}) ---')
        for line in signal[:40]:
            print(line[:300])
    if args.tail > 0:
        print(f'--- last {min(args.tail, len(lines))} lines ---')
        for line in lines[-args.tail:]:
            print(line[:300])
    sys.exit(proc.returncode)


if __name__ == '__main__':
    main()
