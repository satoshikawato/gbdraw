"""Compare raw Web draft and canonical scalar admission without rendering."""
import copy
import json
import platform
import subprocess
import tempfile
from pathlib import Path
from gbdraw.session import materialize_session
from gbdraw.session_request_codec import decode_canonical_request

root = Path(__file__).resolve().parent
repo = root.parents[3]
seed = json.loads((repo / 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json').read_text())
cases = json.loads((root / 'scalar-fixtures.json').read_text())
rows = []


def decode(value):
    if isinstance(value, dict):
        return float(value['$number']) if '$number' in value else {key: decode(item) for key, item in value.items()}
    return value


with tempfile.TemporaryDirectory(prefix='gbdraw-619-s00-native-') as directory:
    with materialize_session(seed, output_directory=directory) as session:
        for sample in cases:
            for field in ('width', 'radius'):
                for representation in ('rawDraft', 'canonical') if sample['valid'] else ('rawDraft',):
                    request = copy.deepcopy(seed['renderRequest'])
                    row = next(r for r in request['diagramOptions']['tracks']['circularTrackSlots'] if r['id'] == 'gc_content')
                    value = sample['input'] if representation == 'rawDraft' else sample['canonical']
                    row[field] = decode(value)
                    try:
                        decode_canonical_request(request, resource_paths=session.resource_paths, output_directory=directory)
                        result = {'accepted': True}
                    except Exception as error:
                        result = {'accepted': False, 'error': str(error)}
                    if representation == 'canonical':
                        assert result['accepted'], (sample['name'], field, result)
                    rows.append({'name': sample['name'], 'field': field, 'representation': representation, 'input': value, 'nativeCanonicalRequest': result})
print(json.dumps({'sourceSha': subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(), 'pythonVersion': platform.python_version(), 'observations': rows, 'limits': ['Native request conversion only; no render, geometry, CLI replay or export.', 'Acceptance of raw Web draft by native decoder does not authorize widening canonical Web request output.']}, indent=2))
