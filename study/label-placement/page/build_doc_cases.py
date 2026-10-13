#!/usr/bin/env python3
"""Add the Circular commands of the current Gallery and tutorials (origin/dev tree in work/dev-tree) to extra_cases.json.

Gallery: examples.json `command`; tutorials: `<!-- executable:T-CLI-nn:start -->` blocks in docs/TUTORIALS/*.md.
Each command becomes one case: `-o` / `-f` removed, input file names replaced by absolute paths in the dev tree.
"""
import json, re, shlex
from pathlib import Path

T = Path('/home/kawato/gbdraw-baselines/label-placement-20261011')
DEV = T / 'work' / 'dev-tree'
CASES = T / 'page' / 'extra_cases.json'


def index(*roots):
    found = {}
    for root in roots:
        for path in sorted(root.rglob('*')):
            if path.is_file():
                found.setdefault(path.name, path)
    return found


def to_args(command, files):
    tokens = shlex.split(command)
    assert tokens[:2] == ['gbdraw', 'circular'], tokens[:2]
    args, skip, out = ['circular'], False, None
    for i, tok in enumerate(tokens[2:], 2):
        if skip:
            skip = False
            continue
        if tok in ('-o', '--output', '-f', '--format'):
            if tok in ('-o', '--output'):
                out = tokens[i + 1]
            skip = True
            continue
        if tok not in files and '/' in tok and not tok.startswith('-'):
            hits = [p for p in files.values() if str(p).endswith('/' + tok)]
            tok = str(hits[0]) if hits else tok
        args.append(str(files[tok]) if tok in files else tok)
    missing = [a for a in args if re.fullmatch(r'[A-Za-z0-9_./-]+\.(gbk|gb|gbff|tsv|out|gff3|fna|fasta|txt|toml|gz)', a) and not a.startswith('/')]
    return args, out, missing


data = json.loads(CASES.read_text())
cases = data['cases']
added = []
gallery_files = index(DEV / 'gbdraw/web/gallery/files', DEV / 'gbdraw/web/tutorial-data')
for ex in json.loads((DEV / 'gbdraw/web/gallery/examples.json').read_text()):
    cmd = ex.get('command', '')
    if not cmd.startswith('gbdraw circular'):
        continue
    args, out, missing = to_args(cmd, gallery_files)
    cid = 'gal_' + re.sub(r'[^a-z0-9]+', '_', ex['id'].lower()).strip('_')
    cases[cid] = {'record': f"Gallery: {ex['id']}", 'title': re.sub('<[^>]+>', '', ex['title']), 'source': 'gallery', 'args': args, 'missing': missing}
    added.append(cid)
tutorial_files = index(DEV / 'gbdraw/web/tutorial-data')
for md in sorted((DEV / 'docs/TUTORIALS').glob('*.md')):
    text = md.read_text()
    heading = re.search(r'^# (.+)$', text, re.M).group(1)
    for sid, block in re.findall(r'<!-- executable:(T-CLI-\d+):start -->\s*```(?:bash|sh)?\n(.*?)```', text, re.S):
        for cmd in re.split(r'\n(?=gbdraw )', block.replace('\\\n', ' ').strip()):
            cmd = ' '.join(cmd.split('\n')[0].split()) if cmd.startswith('gbdraw') else ''
            if not cmd.startswith('gbdraw circular'):
                continue
            args, out, missing = to_args(cmd, tutorial_files)
            cid = 'doc_' + sid.lower().replace('-', '') + '_' + re.sub(r'[^a-z0-9]+', '_', (out or 'out').lower()).strip('_')
            cases[cid] = {'record': f'Tutorial {sid}: {md.name} -> {out}', 'title': heading, 'source': 'tutorial', 'args': args, 'missing': missing}
            added.append(cid)
# tables that a tutorial tells the reader to save by hand: written by the scenario runner's own table source
import sys
sys.path.insert(0, str(DEV / 'docs' / 'recipes'))
import run_cli_scenarios
for cid in added:
    case = cases[cid]
    if not any(a.startswith('tables/') for a in case['args']):
        continue
    sid = re.search(r'T-CLI-\d+', case['record']).group(0)
    workdir = T / 'work' / 'doc-work' / sid
    workdir.mkdir(parents=True, exist_ok=True)
    run_cli_scenarios._materialize_generated_tables(sid, workdir)
    case['args'] = [str(workdir / a) if a.startswith('tables/') else a for a in case['args']]
    case['missing'] = [a for a in case['missing'] if not (workdir / a).exists()]
for cid in added:
    args, drop = [], False
    for a in cases[cid]['args']:
        if drop:
            drop = False
            continue
        if a == '--session_output':
            drop = True
            cases[cid]['missing'] = [m for m in cases[cid]['missing'] if not m.endswith('.gbdraw-session.json.gz')]
            continue
        args.append(a)
    cases[cid]['args'] = args
CASES.write_text(json.dumps(data, ensure_ascii=False, indent=1))
for cid in added:
    c = cases[cid]
    opts = [a for a in c['args'] if a.startswith('-')]
    print(cid, '|', c['title'][:60], '| labels:', [c['args'][i + 1] if i + 1 < len(c['args']) and not c['args'][i + 1].startswith('-') else 'flag' for i, a in enumerate(c['args']) if a == '--labels'],
          '| placement:', [c['args'][i + 1] for i, a in enumerate(c['args']) if a == '--label_placement'], '| missing:', c['missing'])
