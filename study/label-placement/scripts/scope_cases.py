#!/usr/bin/env python3
"""Write scripts/scope_cases.json: the case list of the scope study (LP-D48 Q10). Run once; edit here, not the JSON."""
import json, pathlib

T = pathlib.Path(__file__).resolve().parent.parent
WT = '/home/kawato/gbdraw-work/.worktrees/lp-scope'
TI = f'{WT}/tests/test_inputs'
EXI = f'{WT}/examples'
EX = '/home/kawato/gbdraw-baselines/gallery-video-20261010/s4/len100k/extra-inputs'
IN = str(T / 'work' / 'inputs')
extra = json.loads((T / 'page' / 'extra_cases.json').read_text())['cases']
std = json.loads((T / 'scripts' / 'cases.json').read_text())
GBK = {'p1': f'{EX}/NC_005856.1.gb', 't5': f'{EX}/NC_005859.1.gb', 'n4': f'{EX}/NC_008720.1.gb', 'epifagus': f'{EX}/NC_001568.1.gb',
       'tobacco': f'{TI}/NC_001879.gbk', 'mito': f'{TI}/HmmtDNA.gbk', 'lambda': f'{TI}/NC_001416.gb', 'mjenmv': f'{TI}/MjeNMV.gbk',
       'wssv': f'{IN}/AP027280.gb', 'plasmid37': f'{IN}/NZ_AP024092.1.gb', 'plasmid153': f'{IN}/NZ_CP125877.1.gb',
       'mg1655': f'{TI}/MG1655.gbk'}
cases = {}


def both(args):
    a = list(args)
    if '--labels' in a:
        i = a.index('--labels')
        if i + 1 < len(a) and not a[i + 1].startswith('-'):
            a[i + 1] = 'both'
        else:
            a.insert(i + 1, 'both')
    else:
        a += ['--labels', 'both']
    return a


def setopt(args, opt, value):
    a = list(args)
    if opt in a:
        a[a.index(opt) + 1] = value
    else:
        a += [opt, value]
    return a


# ---- Circular, horizontal, --labels both ----
COMMON = ['--separate_strands', '--track_type', 'tuckin', '--labels', 'both', '--legend', 'right']
cases['c_mjenmv_ref_both'] = ['circular', '--gbk', GBK['mjenmv'], '--separate_strands', '--labels', 'both', '--legend', 'none']
cases['c_mjenmv_ex_both'] = both(extra['ex_mjenmv_labels']['args'])
cases['c_wssv_ex_both'] = both(extra['ex_wssv_labels']['args'])
for key, font in [('p1', None), ('p1', 14), ('t5', None), ('t5', 14), ('n4', None), ('tobacco', None), ('tobacco', 14),
                  ('mito', None), ('mito', 20), ('plasmid37', None), ('plasmid153', None), ('plasmid153', 14),
                  ('lambda', None), ('epifagus', None), ('mjenmv', None), ('wssv', None)]:
    name = f'c_{key}_both' + (f'_f{font}' if font else '_def')
    cases[name] = ['circular', '--gbk', GBK[key], *COMMON] + (['--label_font_size', str(font)] if font else [])
cases['c_pbr322_both_def'] = both(extra['pbr322_f14']['args'])
cases['c_mg1655rna_both_def'] = both(extra['mg1655_rna_f8']['args'])
cases['c_lambda_mid_both_def'] = ['circular', '--gbk', GBK['lambda'], '--track_type', 'middle', '--labels', 'both']   # docs/RECIPES.md
cases['c_tobacco_mid_both_def'] = ['circular', '--gbk', GBK['tobacco'], '--track_type', 'middle', '--labels', 'both']
cases['c_doc_tcli03_both'] = list(extra['doc_tcli03_mitochondrial_features_highlighted']['args'])
cases['c_doc_tcli06_hz_both'] = list(extra['doc_tcli06_cli_annotated_chloroplast_hz']['args'])
cases['c_gal_tobacco_hz_both'] = list(extra['gal_tobacco_chloroplast_hz']['args'])
cases['c_gal_tobacco_hz_in1_both'] = list(extra['gal_tobacco_chloroplast_hz_in1']['args'])
cases['c_gal_atskew_both'] = both(extra['gal_hmmtdna_atskew']['args'])

# ---- Circular, radial ----
cases['r_mjenmv_ref'] = ['circular', '--gbk', GBK['mjenmv'], '--labels', 'both', '--label_placement', 'radial',
                         '--label_rendering', 'external_only', '--legend', 'none']
cases['r_gal_tobacco'] = list(extra['gal_tobacco_chloroplast']['args'])
cases['r_doc_tcli06'] = list(extra['doc_tcli06_cli_annotated_chloroplast']['args'])
for key, font in [('p1', 14), ('t5', None), ('tobacco', 14), ('n4', None), ('wssv', None)]:
    name = f'r_{key}_both' + (f'_f{font}' if font else '_def')
    cases[name] = ['circular', '--gbk', GBK[key], *COMMON, '--label_placement', 'radial'] + (['--label_font_size', str(font)] if font else [])

# ---- Linear ----
GAL = f'{WT}/gbdraw/web/gallery/sessions'
cases['l_gal_lambda'] = ['linear', '--session', f'{GAL}/lambda_basic_linear.gbdraw-session.json']
cases['l_gal_bgc'] = ['linear', '--session', f'{GAL}/BGC0000708-BGC0000713.gbdraw-session.json']
cases['l_mjenmv_ref'] = ['linear', '--gbk', GBK['mjenmv'], '--show_labels', 'all', '--legend', 'none']
cases['l_mjenmv_ref_ss'] = ['linear', '--gbk', GBK['mjenmv'], '--separate_strands', '--show_labels', 'all', '--legend', 'none']
cases['l_mjenmv_multi'] = ['linear', '--gbk', GBK['mjenmv'], f'{TI}/MelaMJNV.gb' if pathlib.Path(f'{TI}/MelaMJNV.gb').exists() else f'{EXI}/MelaMJNV.gb',
                           '--show_labels', 'all', '--legend', 'none']
cases['l_mjenmv_blast'] = cases['l_mjenmv_multi'] + ['-b', f'{EXI}/MjeNMV.MelaMJNV.tblastx.out', '--align_center']
chain = ['MjeNMV', 'MelaMJNV', 'PemoMJNVA', 'PeseMJNV']
cases['l_majani4_blast'] = ['linear', '--gbk', *[f'{EXI}/{c}.gb' for c in chain], '-b',
                            *[f'{EXI}/{a}.{b}.tblastx.out' for a, b in zip(chain, chain[1:])], '--align_center', '--show_labels', 'all']
cases['l_majani4_blast_f8'] = cases['l_majani4_blast'] + ['--label_font_size', '8']
cases['l_majani4_blast_ss'] = cases['l_majani4_blast'] + ['--separate_strands']
hep = ['AP027078', 'AP027131', 'AP027133', 'AP027132', 'NZ_CP006932']
cases['l_hepato5_blast'] = ['linear', '--gbk', *[f'{TI}/{c}.gb' for c in hep], '-b',
                            *[f'{TI}/{a}_{b}.tblastx.out' for a, b in zip(hep, hep[1:])], '--align_center', '--separate_strands', '--show_labels', 'all']
cases['l_hepato5_above'] = cases['l_hepato5_blast'] + ['--label_placement', 'above_feature']
for key, fonts in [('t5', [None, 8, 10]), ('tobacco', [None, 8, 10]), ('wssv', [None, 8, 10]), ('p1', [None, 32]),
                   ('lambda', [None, 32]), ('mjenmv', [8]), ('mg1655', [None])]:
    for font in fonts:
        name = f'l_{key}' + (f'_f{font}' if font else '_def')
        cases[name] = ['linear', '--gbk', GBK[key], '--show_labels', 'all'] + (['--label_font_size', str(font)] if font else [])
    cases[f'l_{key}_ss'] = ['linear', '--gbk', GBK[key], '--separate_strands', '--show_labels', 'all']
    cases[f'l_{key}_above'] = ['linear', '--gbk', GBK[key], '--separate_strands', '--show_labels', 'all', '--label_placement', 'above_feature']
for key in ('t5', 'p1', 'lambda', 'tobacco'):
    cases[f'l_{key}_above_r45'] = cases[f'l_{key}_above'] + ['--label_rotation', '45']
for key in ('t5', 'wssv', 'tobacco'):
    cases[f'l_{key}_resolve'] = ['linear', '--gbk', GBK[key], '--show_labels', 'all', '--resolve_overlaps']
cases['l_lambda_below'] = ['linear', '--gbk', GBK['lambda'], '--show_labels', 'all', '--track_layout', 'below']
cases['l_t5_extonly'] = ['linear', '--gbk', GBK['t5'], '--show_labels', 'all', '--label_rendering', 'external_only']

(T / 'scripts' / 'scope_cases.json').write_text(json.dumps({'_doc': __doc__, 'cases': cases}, indent=1))
print(len(cases), 'cases')
for k, v in cases.items():
    for a in v:
        if a.startswith('/') and not pathlib.Path(a).exists():
            print('MISSING', k, a)
