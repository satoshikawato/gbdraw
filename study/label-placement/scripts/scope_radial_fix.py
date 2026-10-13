#!/usr/bin/env python3
"""Measurement only: run scope_run.py with circular_radial._place_side compiled from its source with the inner
preflight threshold 1e-6 raised to 1e-3 (the threshold at which assemble stops growing the radius). usage: as scope_run.py"""
import inspect, runpy, sys, pathlib, textwrap
import gbdraw.labels.circular_radial as cr
src = textwrap.dedent(inspect.getsource(cr._place_side))
assert 'required_growth > 1e-6' in src
exec(compile(src.replace('required_growth > 1e-6', 'required_growth > 1e-3'), cr.__file__, 'exec'), cr.__dict__)
sys.argv = [str(pathlib.Path(__file__).with_name('scope_run.py'))] + sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
