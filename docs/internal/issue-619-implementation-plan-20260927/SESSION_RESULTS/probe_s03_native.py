"""Replay actual S03 downloads through the existing native CLI; never edit inputs."""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import math
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

GEOMETRY_TAGS = {
    'svg', 'g', 'path', 'circle', 'ellipse', 'rect', 'line', 'polyline',
    'polygon', 'text', 'textPath',
}
GEOMETRY_ATTRIBUTES = {
    'viewBox', 'transform', 'd', 'cx', 'cy', 'r', 'rx', 'ry', 'x', 'y',
    'x1', 'x2', 'y1', 'y2', 'points', 'width', 'height', 'font-size', 'startOffset',
}
NUMBER = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?')
ABS_TOLERANCE = 1e-10


def read_document(path: Path) -> dict:
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data)


def svg_geometry(content: str) -> list:
    result = []
    for node in ET.fromstring(content).iter():
        name = node.tag.rsplit('}', 1)[-1]
        if name in GEOMETRY_TAGS:
            result.append((name, {
                k: v for k, v in node.attrib.items() if k in GEOMETRY_ATTRIBUTES
            }, ''.join(node.itertext()) if name in {'text', 'textPath'} else None))
    return result


def compare_svg(expected: str, actual: str) -> dict:
    left, right = svg_geometry(expected), svg_geometry(actual)
    assert len(left) == len(right), (len(left), len(right))
    numeric_differences = []
    for index, ((tag_a, attrs_a, text_a), (tag_b, attrs_b, text_b)) in enumerate(zip(left, right)):
        assert (tag_a, text_a) == (tag_b, text_b), (index, tag_a, tag_b, text_a, text_b)
        assert attrs_a.keys() == attrs_b.keys(), (index, attrs_a, attrs_b)
        for key in attrs_a:
            a, b = attrs_a[key], attrs_b[key]
            if a == b:
                continue
            assert NUMBER.sub('#', a) == NUMBER.sub('#', b), (index, key, a, b)
            x, y = NUMBER.findall(a), NUMBER.findall(b)
            assert len(x) == len(y)
            for u, v in zip(x, y):
                assert math.isclose(float(u), float(v), rel_tol=0, abs_tol=ABS_TOLERANCE), (index, key, u, v)
                numeric_differences.append(abs(float(u) - float(v)))
    return {
        'elementsCompared': len(left), 'textAndTransformStructureEqual': True,
        'absoluteTolerancePx': ABS_TOLERANCE,
        'maxNumericDifference': max(numeric_differences, default=0),
        'differentNumericTokens': len(numeric_differences),
        'byteIdentity': expected.encode() == actual.encode(),
    }


def scientific_geometry(value: dict) -> dict:
    data = json.loads(json.dumps(value))
    for record in data['records']:
        record.pop('resultName', None)  # Explicit CLI output override, not geometry.
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    source = args.session.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    document = read_document(source)
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    mode = document['renderRequest']['mode']
    prefix = output / 'native'
    sidecar = output / 'native.gbdraw-session.json'
    command = [sys.executable, '-m', 'gbdraw.cli', mode, '--session', str(source),
               '--output', str(prefix), '--format', 'svg', '--save_session',
               '--session_output', str(sidecar)]
    # A clean output directory and the installed/explicit PYTHONPATH package are
    # the only context. Embedded biological resources come from the real download.
    run = subprocess.run(command, cwd=output, capture_output=True, text=True, check=False)
    (output / 'stdout.log').write_text(run.stdout)
    (output / 'stderr.log').write_text(run.stderr)
    assert run.returncode == 0, run.stderr
    native = read_document(sidecar)
    assert native['renderRequest']['diagramOptions']['tracks'] == document['renderRequest']['diagramOptions']['tracks']
    for expected, actual in zip(document['renderRequest']['records'], native['renderRequest']['records'], strict=True):
        original = dict(expected)
        replayed = dict(actual)
        selector = original.pop('selector')
        replay_selector = replayed.pop('selector')
        assert original == replayed
        resource_id = expected['source']['resourceId']
        assert base64.b64decode(document['resources'][resource_id]['data']) == base64.b64decode(native['resources'][resource_id]['data'])
        if selector != replay_selector:
            # Native resolves a selector over a one-record GenBank source and
            # omits that redundant selector in its sidecar. Prove the identity.
            from Bio import SeqIO
            import io
            records = list(SeqIO.parse(io.StringIO(base64.b64decode(document['resources'][resource_id]['data']).decode()), 'genbank'))
            assert replay_selector is None and selector['kind'] == 'recordId'
            assert [record.id for record in records] == [selector['value']]
    assert scientific_geometry(native['runMetadata']['trackSlotGeometry']) == scientific_geometry(document['runMetadata']['trackSlotGeometry'])
    svg = prefix.with_suffix('.svg')
    comparison = compare_svg(document['results'][0]['content'], svg.read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest() == source_hash
    result = {
        'session': str(source), 'sessionSha256': source_hash, 'command': command,
        'cwd': str(output), 'exit': run.returncode, 'svg': str(svg),
        'svgSha256': hashlib.sha256(svg.read_bytes()).hexdigest(), 'sidecar': str(sidecar),
        'canonicalValuesAndUnitsExact': True, 'resolvedGeometryExactExceptOutputName': True,
        'draftUsedForReplay': False, 'svgComparison': comparison,
        'nativeEncoding': 'CLI sidecar re-encodes file-backed color/qualifier tables as canonical TSV; default-valued options and redundant single-record selectors may be omitted. Biological bytes and track scalars remain exact.',
    }
    (output / 'comparison.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
