"""Exercise the exact embedded Python render/helper boundary for Node transport tests."""
from __future__ import annotations

import base64
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.api import CircularDiagramOptions, CircularDiagramRequest, ColorOptions, InMemoryRecordSource, RecordInput
from gbdraw.session import build_session_document


def python_helpers():
    source = (Path(__file__).resolve().parents[3] / 'gbdraw/web/js/app/python-helpers.js').read_text()
    namespace = {}
    exec(source.split('`', 1)[1].rsplit('`', 1)[0], namespace)
    return namespace


def render_failure(namespace, pattern):
    record = SeqRecord(Seq('ATGC' * 25), id='PRIVATE_RECORD_SENTINEL')
    record.annotations['molecule_type'] = 'DNA'
    document = build_session_document(CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        options=CircularDiagramOptions(colors=ColorOptions(color_table=DataFrame([
            dict(feature_type='CDS', qualifier_key='product', value='VALID_PLACEHOLDER', color='#ff0000', caption='')
        ]))),
    )).to_dict()
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
            path.write_bytes(data)
            paths[resource_id] = str(path)
        result = namespace['run_canonical_request_wrapper'](
            json.dumps(document['renderRequest']), json.dumps(paths), str(workspace))
        assert not workspace.exists()
        return result


def main():
    payload = json.load(sys.stdin)
    namespace = python_helpers()
    if payload.get('render'):
        result = render_failure(namespace, payload['pattern'])
    else:
        result = json.loads(namespace['call_web_json_helper'](payload['helper'], *payload['args']))
    print(json.dumps(result))


if __name__ == '__main__':
    main()
