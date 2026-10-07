"""Exercise the exact embedded Python render/helper boundary for Node transport tests."""
from __future__ import annotations

import base64
import builtins
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

# Resolve native producers from the same checkout as the embedded Web helpers.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from gbdraw.api import CircularDiagramOptions, CircularDiagramRequest, ColorOptions, InMemoryRecordSource, RecordInput
from gbdraw import exceptions
from gbdraw.session import build_session_document
from gbdraw.web_support import request_render

# A native failure inside a real engine phase: the phase function it replaces.
ENGINE_PHASES = {'render': 'render_request', 'result-admission': 'build_feature_catalog_item'}


def python_helpers():
    source = (Path(__file__).resolve().parents[3] / 'gbdraw/web/js/app/python-helpers.js').read_text()
    namespace = {}
    exec(source.split('`', 1)[1].rsplit('`', 1)[0], namespace)
    return namespace


def render_failure(namespace, pattern, display_start=None, genbank_text=None, config_overrides=None):
    record = SeqRecord(Seq('ATGC' * 25), id='PRIVATE_RECORD_SENTINEL')
    record.annotations['molecule_type'] = 'DNA'
    document = build_session_document(CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=CircularDiagramOptions(colors=ColorOptions(color_table=DataFrame([
            dict(feature_type='CDS', qualifier_key='product', value='VALID_PLACEHOLDER', color='#ff0000', caption='')
        ]))),
    )).to_dict()
    if display_start is not None:
        document['renderRequest']['records'][0]['display']['startCoordinate'] = display_start
    if config_overrides is not None:
        document['renderRequest']['diagramOptions']['configOverrides'].update(config_overrides)
    with TemporaryDirectory(prefix='structured-error-') as root:
        workspace = Path(root) / 'gbdraw-web-render-1'
        resources = workspace / 'resources'
        resources.mkdir(parents=True)
        (workspace / '.gbdraw-worker-render-workspace').touch()
        paths = {}
        for index, (resource_id, entry) in enumerate(document['resources'].items()):
            path = resources / f'{index}.bin'
            data = base64.b64decode(entry['data'])
            if resource_id == 'colors-color-table':
                data = data.replace(b'VALID_PLACEHOLDER', pattern.encode('utf-8'))
            if genbank_text is not None and resource_id.endswith('-genbank'):
                data = genbank_text.encode('utf-8')
            path.write_bytes(data)
            paths[resource_id] = str(path)
        result = namespace['run_canonical_request_wrapper'](
            json.dumps(document['renderRequest']), json.dumps(paths), str(workspace))
        assert not workspace.exists()
        return result


def engine_phase_failure(namespace, name, stage):
    exception = getattr(exceptions, name, None) or getattr(builtins, name)

    def fail(*_args, **_kwargs):
        raise exception('PRIVATE_EXCEPTION_SENTINEL')

    phase = ENGINE_PHASES[stage]
    original = getattr(request_render, phase)
    setattr(request_render, phase, fail)
    try:
        return render_failure(namespace, 'VALID_PLACEHOLDER')
    finally:
        setattr(request_render, phase, original)


def main():
    payload = json.load(sys.stdin)
    namespace = python_helpers()
    if payload.get('raise'):
        result = [engine_phase_failure(namespace, name, stage) for name, stage in payload['raise']]
    elif any(payload.get(key) is not None for key in ('displayStart', 'genbankText', 'configOverrides')):
        result = render_failure(namespace, 'VALID_PLACEHOLDER', payload.get('displayStart'), payload.get('genbankText'),
                                payload.get('configOverrides'))
    elif payload.get('render'):
        result = render_failure(namespace, payload['pattern'])
    else:
        result = json.loads(namespace['call_web_json_helper'](payload['helper'], *payload['args']))
    print(json.dumps(result))


if __name__ == '__main__':
    main()
