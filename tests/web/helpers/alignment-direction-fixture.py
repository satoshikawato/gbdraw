"""Build the S01 source fixture as a real browser Session, without LOSAT."""
from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult  # noqa: E402
from gbdraw.api import build_session_document, materialize_session, save_session_document  # noqa: E402
from gbdraw.web_support.request_render import render_canonical_web_request  # noqa: E402
from tests.test_alignment_direction_projection import build_fixture  # noqa: E402


def main(directory: Path) -> None:
    request, helper, _context, _paths, _response, _sources = build_fixture(directory)
    # The negative-strand reference is the only left-facing selected anchor.
    # UI fixture uses full records; native S01 tests retain the crop geometry cases.
    records = tuple(replace(record, region=None, presentation=replace(record.presentation, label=record.record_key, reverse_complement=False),
                            record_key=f'record-{index + 1}')
                    for index, record in enumerate(request.records) if record.record_key != 'unusable')
    members = []
    for member in helper['members']:
        if member['anchor']['recordKey'] == 'unusable':
            continue
        key = member['anchor']['recordKey']
        members.append(OrthogroupMember(
            orthogroup_id='group', protein_id=key + '-protein',
            record_index=next(i for i, record in enumerate(request.records) if record.record_key == key),
            feature_index=0, record_id=key, label=key + '-protein',
            start=member['sourceStart'], end=member['sourceEnd'], strand=member['sourceStrand'],
            feature_svg_id=None, source_protein_id=key + '-protein'))
    request = replace(request, records=records, layout=replace(request.layout,
        record_translations=tuple(replace(entry, record_key=f'record-{index + 1}')
            for index, entry in enumerate(request.layout.record_translations) if entry.record_key != 'unusable')),
        options=replace(request.options,
        orthogroups=OrthogroupResult({'group': members}, {m.protein_id: m for m in members})))
    document = build_session_document(request)
    with materialize_session(document, output_directory=directory / 'decoded') as materialized:
        web = render_canonical_web_request(document.to_dict()['renderRequest'],
            resource_paths=materialized.resource_paths, output_directory=directory / 'web')
    save_session_document(directory / 'directions.gbdraw-session.json', request, adjunct={
        'results': web['results'], 'editorState': {'featureCatalog': web['metadata']['featureCatalog']}})


if __name__ == '__main__':
    main(Path(sys.argv[1]))
