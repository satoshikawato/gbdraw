"""case_args / apply_variant of fit_sweep.py, importable (generated from it)."""

import argparse, json, os, pathlib, subprocess, time
from concurrent.futures import ThreadPoolExecutor

HERE = pathlib.Path(__file__).resolve().parent
TASK = HERE.parent
ROOTS = {'EX': '/home/kawato/gbdraw-baselines/gallery-video-20261010/s4/len100k/extra-inputs',
         'IN': str(TASK / 'work' / 'inputs')}
ROOTS['TI'] = '/home/kawato/gbdraw-work/.worktrees/lp-fit/tests/test_inputs'
extra = json.loads((TASK / 'page' / 'extra_cases.json').read_text())['cases']
std = json.loads((HERE / 'cases.json').read_text())


def case_args(case):
    if case in extra:
        return list(extra[case]['args'])
    c = std['cases'][case]
    root, rel = c['gbk'].split('/', 1)
    return ['circular', '--gbk', f'{ROOTS[root]}/{rel}', *std['common'], '--label_font_size', str(c['font'])]


def set_opt(args, opt, value):
    if opt in args:
        i = args.index(opt); args[i + 1] = value
    else:
        args += [opt, value]


def apply_variant(args, variant):
    for part in variant.split('+'):
        if part == 'cur':
            continue
        if part.startswith('xy'):
            set_opt(args, '--outer_label_x_radius_offset', part[2:]); set_opt(args, '--outer_label_y_radius_offset', part[2:])
        elif part[0] == 'x':
            set_opt(args, '--outer_label_x_radius_offset', part[1:])
        elif part[0] == 'y':
            set_opt(args, '--outer_label_y_radius_offset', part[1:])
        elif part[0] == 'f':
            set_opt(args, '--label_font_size', part[1:])
        elif part[0] == 's':
            set_opt(args, '--circular_label_spacing', part[1:])
        else:
            raise SystemExit(f'unknown variant part {part!r}')
    return args


