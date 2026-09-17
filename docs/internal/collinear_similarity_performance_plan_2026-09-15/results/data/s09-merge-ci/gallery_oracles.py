"""Compare historical Gallery science on an explicit source tree and interpreter.

This runs saved raw analysis only, with no timing or new searches. The output
contains exact digests; optional canonical files are local diagnostic artifacts.
"""
import argparse
from dataclasses import replace
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--canonical-dir', type=Path)
    args = parser.parse_args()
    root = args.source_root.resolve()
    runner_path = Path(__file__).resolve().parents[6] / 'tools/benchmark_protein_comparison.py'
    spec = importlib.util.spec_from_file_location('gallery_oracle_runner', runner_path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    sys.path.insert(0, str(root))
    pc, cc = runner.load_source(root)
    import pandas as pd

    source_hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted((root / 'gbdraw').rglob('*.py'))}
    report = {
        'python': platform.python_version(),
        'dependencies': {n: importlib.metadata.version(n) for n in ('pandas', 'numpy', 'biopython')},
        'sourcePythonSha256': hashlib.sha256(json.dumps(source_hashes, sort_keys=True).encode()).hexdigest(),
        'runnerSha256': hashlib.sha256(runner_path.read_bytes()).hexdigest(),
        'pandasInferString': False,
        'cases': {},
    }
    with pd.option_context('future.infer_string', False):
        for name in ('gallery-collinear', 'gallery-orthogroup'):
            inventory, stages = runner.build_case(root, pc, cc, name, runner.SEED)
            result = stages['post_search']()
            if hasattr(pc, 'materialize_ortholog_paths'):
                result = replace(result, orthogroups=pc.materialize_ortholog_paths(result.orthogroups))
            canonical = runner.json_bytes(runner.canonical(result))
            report['cases'][name] = {
                'semanticSha256': runner.digest(canonical),
                'canonicalBytes': len(canonical),
                'sessionSha256': inventory['sessionSha256'],
            }
            if args.canonical_dir:
                args.canonical_dir.mkdir(parents=True, exist_ok=True)
                (args.canonical_dir / f'{name}.json.gz').write_bytes(gzip.compress(canonical, mtime=0))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
