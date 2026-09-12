"""Portable settings-only document controls; the fixture is a real Web download."""
import copy
import gzip
import json
from pathlib import Path

import pytest

from gbdraw.api import load_session_document, materialize_session, session_to_request
from gbdraw.session import SessionConversionError, SessionError


@pytest.fixture
def settings():
    path = Path(__file__).parent / 'fixtures/sessions/settings-only.v42.json.gz'
    return json.loads(gzip.decompress(path.read_bytes()))


def test_settings_only_reader_and_nonrenderable_edge(settings, tmp_path):
    document = load_session_document(settings)
    assert document.version == 42
    assert not document.has_canonical_request
    assert document.to_dict()['config'] == settings['config']
    assert settings['config']['form']['labels_mode'] == 'both'
    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionConversionError, match='Settings-only Session has no biological render request'):
            session_to_request(materialized)


def test_settings_only_auxiliary_resource_materializes(settings, tmp_path):
    import base64
    content = b'CDS\tgene,product,locus_tag\r\n'
    settings['resources']['priority'] = {
        'kind': 'web-file', 'name': 'priority.tsv', 'type': 'text/tab-separated-values',
        'size': len(content), 'lastModified': 123, 'encoding': 'base64',
        'data': base64.b64encode(content).decode(),
    }
    settings['webFiles']['bindings']['qualifier_priority'] = {
        'resourceId': 'priority', 'name': 'priority.tsv', 'type': 'text/tab-separated-values', 'lastModified': 123,
    }
    document = load_session_document(settings)
    with materialize_session(document, output_directory=tmp_path) as materialized:
        assert materialized.resource_paths['priority'].read_bytes() == content
        assert document.to_dict()['webFiles'] == settings['webFiles']


@pytest.mark.parametrize('damage,reason', [
    ('old_version', 'requires a canonical renderRequest'),
    ('future_version', 'newer than'),
    ('missing_request', 'requires a canonical renderRequest'),
    ('config', 'active Web configuration'),
    ('binding', 'missing resource'),
    ('binding_schema', 'binding schema'),
    ('source', 'biological sources'),
    ('result', 'editorState.featureCatalog'),
    ('provenance', 'committed render artifacts'),
])
def test_settings_only_rejects_inconsistent_documents(settings, damage, reason):
    document = copy.deepcopy(settings)
    if damage == 'old_version':
        document['version'] = 41
    elif damage == 'future_version':
        document['version'] = 999
    elif damage == 'missing_request':
        del document['renderRequest']
    elif damage == 'config':
        document['config']['form'] = []
    elif damage == 'binding':
        document['webFiles']['bindings']['whitelist'] = {
            'resourceId': 'absent', 'name': 'list.tsv', 'type': '', 'lastModified': 0,
        }
    elif damage == 'binding_schema':
        document['webFiles']['bindings']['schema'] = 99
    elif damage == 'source':
        full = load_session_document(Path(__file__).parent / 'fixtures/sessions/single.v41-bindings1.json').to_dict()
        document['resources'] = full['resources']
        document['webFiles'] = full['webFiles']
    elif damage == 'result':
        document['results'] = [{'name': 'old.svg', 'content': '<svg/>'}]
    elif damage == 'provenance':
        document['cliInvocation'] = {}
    with pytest.raises(SessionError, match=reason):
        load_session_document(document)


def test_cli_rejects_settings_only_without_creating_output(settings, tmp_path):
    from gbdraw.circular import circular_main
    from gbdraw.exceptions import ValidationError
    source = tmp_path / 'settings.json'
    source.write_text(json.dumps(settings))
    with pytest.raises(ValidationError, match='Settings-only Session has no biological render request'):
        circular_main(['--session', str(source), '-o', str(tmp_path / 'out')])
    assert not (tmp_path / 'out.svg').exists()
