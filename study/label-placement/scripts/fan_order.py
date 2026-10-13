#!/usr/bin/env python3
"""LP-D19: clusters whose straight leaders are not one unbroken run.

The Owner: in a cluster the outermost leaders may bend and the middle ones are straight; straight outside, bent
inside, straight again is not good. Cluster = the definition of work/J/lop.py (touching boxes of one column and one
half, at least one leader turning 3 degrees or more). A leader is straight under `--straight` degrees of turn
(default 1.0). A cluster is "mixed" when a bent leader lies between two straight leaders of the same cluster.

usage: fan_order.py [--straight DEG] [-v] NAME=DIR:SUFFIX ...     (files DIR/<case>__<SUFFIX>.labels.json)
"""
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'work', 'J'))
import lop  # noqa: E402


def mixed(turns, straight=1.0):
    """True when a bent leader lies between two straight ones (turns in stack order)."""
    flat = [abs(t) < straight for t in turns]
    if True not in flat:
        return False
    first, last = flat.index(True), len(flat) - 1 - flat[::-1].index(True)
    return not all(flat[first:last + 1])


def count(path, straight=1.0):
    found = lop.clusters(path)
    bad = [c for c in found if mixed(c['s'], straight)]
    return len(found), bad


if __name__ == '__main__':
    args = sys.argv[1:]
    straight, verbose = 1.0, '-v' in args
    if '--straight' in args:
        straight = float(args[args.index('--straight') + 1]); del args[args.index('--straight'):args.index('--straight') + 2]
    for spec in (a for a in args if '=' in a):
        name, rest = spec.split('=', 1); folder, suffix = rest.rsplit(':', 1)
        total = bad_total = 0
        for path in sorted(glob.glob(os.path.join(folder, f'*__{suffix}.labels.json'))):
            n, bad = count(path, straight)
            total += n; bad_total += len(bad)
            if verbose and bad:
                for c in bad:
                    print(f"  {os.path.basename(path).split('__')[0]}: {c['s']} {c['names'][:6]}")
        print(f'{name}: {bad_total} mixed of {total} clusters')
