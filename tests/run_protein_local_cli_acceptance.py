"""Run with the dedicated installed interpreter outside the checkout, no PYTHONPATH."""
import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import zipfile

import gbdraw
from gbdraw.analysis import collinearity as cc, protein_colinearity as pc
from gbdraw.api import load_session_document, materialize_session, session_to_request, save_session_document, render_request

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert 'PYTHONPATH' not in os.environ
    assert not Path.cwd().is_relative_to(ROOT)
    installed = Path(gbdraw.__file__).resolve().parent
    assert installed.is_relative_to(ROOT / '.venv')
    wheel, = (ROOT / 'gbdraw/web').glob('gbdraw-*-py3-none-any.whl')
    with zipfile.ZipFile(wheel) as archive:
        names = [n for n in archive.namelist() if n.startswith('gbdraw/') and n.endswith('.py')]
        for name in names:
            assert archive.read(name) == (installed.parent / name).read_bytes(), name
    runner = runpy.run_path(str(ROOT / 'tools/benchmark_protein_comparison.py'))
    frozen = runpy.run_path(str(ROOT / 'tests/prototypes/protein_s078.py'))
    artifacts = ROOT / '.venv' / (args.output.stem + '-artifacts')
    artifacts.mkdir(parents=True, exist_ok=True)
    report = {'installedPackage': str(installed), 'packagedPythonFilesMatched': len(names),
              'cwdOutsideCheckout': True, 'pythonpathUnset': True, 'preparedRawGeneration': [], 'replays': []}
    replay_inputs = []
    for name in ('gallery-collinear', 'gallery-orthogroup'):
        _, stages = runner['build_case'](ROOT, pc, cc, name, runner['SEED'])
        result = stages['post_search']()
        with frozen['frozen_collinear_callers']():
            # Similarity's stage calls the public selector on pc directly.
            from unittest.mock import patch
            with patch.object(pc, 'select_rbh_orthogroup_edges_from_directional_hits',
                              frozen['frozen'].select_rbh_orthogroup_edges_from_directional_hits):
                expected = stages['post_search']()
        assert runner['canonical'](result) == runner['canonical'](expected)
        gallery = runner['GALLERIES'][0 if name == 'gallery-collinear' else 1]
        document = load_session_document(ROOT / f'gbdraw/web/gallery/sessions/{gallery}.gbdraw-session.json.gz')
        with tempfile.TemporaryDirectory(prefix='gbdraw-s078-fresh-') as directory:
            with materialize_session(document, output_directory=directory) as materialized:
                request = session_to_request(materialized)
                options = replace(request.options, orthogroups=result.orthogroups,
                                  collinearity_blocks=result if name == 'gallery-collinear' else ())
                request = replace(request, options=options)
                saved = artifacts / f'{name}.gbdraw-session.json'
                save_session_document(saved, request)
                rendered = render_request(request)
                assert rendered.output_paths and all(p.is_file() for p in rendered.output_paths)
                report['preparedRawGeneration'].append({'case': name, 'exactFrozenParity': True,
                    'typedSessionVersion': load_session_document(saved).version,
                    'resultSha256': runner['digest'](runner['json_bytes'](runner['canonical'](result))),
                    'svgSha256': [hashlib.sha256(p.read_bytes()).hexdigest() for p in rendered.output_paths]})
                replay_inputs.append((name, saved))
    replay_inputs.append(('released-v2', ROOT / 'tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz'))
    for label, source in replay_inputs:
        with tempfile.TemporaryDirectory(prefix='gbdraw-s078-replay-') as directory:
            argv = [str(Path(sys.executable).parent / 'gbdraw'), 'linear', '--session', str(source),
                    '-o', 'replayed', '-f', 'svg']
            completed = subprocess.run(argv, cwd=directory, capture_output=True, text=True, check=False)
            log = args.output.with_name(f'{args.output.stem}-{label}.log')
            log.write_text(completed.stdout + completed.stderr)
            assert completed.returncode == 0, log.read_text()
            output = Path(directory) / 'replayed.svg'
            assert output.is_file()
            shutil.copyfile(output, artifacts / f'{label}-replayed.svg')
            report['replays'].append({'case': label, 'sessionVersion': load_session_document(source).version,
                'inputSha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'argv': argv,
                'exitCode': completed.returncode, 'svgSha256': hashlib.sha256(output.read_bytes()).hexdigest()})
            print(label, 'replayed', flush=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
