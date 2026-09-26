"""Issue #598 S01: source-bound direction facts through the real renderer."""
from __future__ import annotations

import copy
import json
import subprocess
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import LinearDiagramOptions, LinearMultiRecordOptions, LinearOutputOptions, LinearRecordTranslation
from gbdraw.api.record_planning import ResolvedRecordCollection, project_similarity_alignment_centers
from gbdraw.api.request_render import plan_linear_request
from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordInput, RecordPresentation
from gbdraw.diagrams.linear.assemble import _final_record_translations
from gbdraw.exceptions import ValidationError
from gbdraw.features.source import build_source_feature_catalog
from gbdraw.io.regions import parse_region_spec
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.similarity_alignment import resolve_similarity_alignment_payload, _anchor

ROOT = Path(__file__).resolve().parents[1]


def project(response, helper_request, intent, binding=None):
    script = """
import {readFileSync} from 'node:fs';
globalThis.window = {Vue:{computed:()=>{},ref:()=>{}}};
const {validateSimilarityAlignmentResolution,projectSimilarityAlignmentDirections} = await import('./gbdraw/web/js/app/similarity-alignment.js');
const [response,request,intent,binding] = JSON.parse(readFileSync(0,'utf8'));
console.log(JSON.stringify(projectSimilarityAlignmentDirections({resolution:validateSimilarityAlignmentResolution(response,request),intent,expectedBinding:binding})));
"""
    result = subprocess.run(["node", "--input-type=module", "-e", script], cwd=ROOT,
                            input=json.dumps([response, helper_request, intent, binding or response["projection"]["binding"]]),
                            capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def build_fixture(tmp_path):
    records = []
    helpers = []
    members = []
    paths = []
    for index, (key, length, strand) in enumerate((('ref', 1000, -1), ('target', 500, 1),
                                                 ('unknown', 800, None), ('skip', 600, -1),
                                                 ('missing', 700, 1), ('unusable', 900, -1))):
        parts = [SimpleLocation(100, 130, strand=strand or 1),
                 SimpleLocation(210, 260, strand=strand if strand is not None else -1)]
        record = SeqRecord(Seq('ACGT' * (length // 4) + 'A' * (length % 4)), id=key, description=key)
        record.annotations.update(molecule_type='DNA', topology='linear')
        record.features = [SeqFeature(CompoundLocation(parts), type='CDS', qualifiers={'protein_id': [key + '-protein']})]
        path = tmp_path / (key + '.gbk')
        SeqIO.write(record, path, 'genbank')
        paths.append(path)
        source = SeqIO.read(path, 'genbank')
        identity = build_source_feature_catalog(source)[0]
        anchor = dict(recordKey=key, biologicalFeatureId=identity.biological_feature_id,
                      sourceFeatureIndex=identity.source_feature_index, stableFeatureSvgId=identity.stable_feature_id)
        region = parse_region_spec('target:51-350:rc') if key == 'target' else (
            parse_region_spec('unusable:401-700') if key == 'unusable' else None)
        presentation = RecordPresentation(reverse_complement=key == 'skip')
        records.append(RecordInput(GenBankInputSource(path), record_key=key, region=region, presentation=presentation))
        helpers.append(dict(recordKey=key, recordLength=length, region=None if region is None else
                            dict(start=region.start, end=region.end, reverseComplement=region.reverse_complement),
                            presentation=dict(reverseComplement=presentation.reverse_complement)))
        if key != 'missing':
            members.append(dict(groupId='group', anchor=anchor, sourceStart=100, sourceEnd=260,
                                sourceStrand=strand, identityIsUnique=True, hidden=False,
                                representative=False, role='member'))
    request = LinearDiagramRequest(records=tuple(records),
        options=LinearDiagramOptions(selected_features_set=('CDS',), output=LinearOutputOptions(legend='none'),
            config_overrides={'canvas.show_gc': False, 'canvas.show_skew': False, 'objects.scale.show': False}),
        layout=LinearMultiRecordOptions(record_translations=tuple(
            LinearRecordTranslation(record.record_key, 37 + index * 11, 9 + index * 7)
            for index, record in enumerate(records))))
    encoded = encode_canonical_request(request)
    paths_by_id = {item.resource_id: str(item.source_path) for item in encoded.resources}
    helper = dict(schema=2, groupId='group', records=helpers, reference=members[0]['anchor'], members=members,
                  directEdges=[], choices=[dict(recordKey='skip', kind='skip', anchor=None)])
    context = dict(canonicalRequest=encoded.payload, orientations=None)
    response = resolve_similarity_alignment_payload(helper, projection=context, resource_paths=paths_by_id)
    return request, helper, context, paths_by_id, response, paths


@pytest.fixture
def fixture(tmp_path):
    return build_fixture(tmp_path)


def _apply(request, projected):
    records = []
    for record, result in zip(request.records, projected['records'], strict=True):
        reverse = result['afterReverseComplement']
        records.append(replace(record, region=replace(record.region, reverse_complement=reverse),
                               presentation=replace(record.presentation, reverse_complement=False)) if record.region else
                       replace(record, presentation=replace(record.presentation, reverse_complement=reverse)))
    translations = tuple(LinearRecordTranslation(row['recordKey'], **row['translation']) for row in projected['records'])
    return replace(request, records=tuple(records), layout=replace(request.layout, record_translations=translations))


@pytest.mark.parametrize('intent,arrows,flips', [
    ({'mode': 'keep'}, [-1, -1, None, None, None, None], []),
    ({'mode': 'right'}, [1, 1, None, None, None, None], ['ref', 'target']),
    ({'mode': 'left'}, [-1, -1, None, None, None, None], []),
    ({'mode': 'custom', 'byRecordKey': {'ref': 'right', 'target': 'keep', 'unknown': 'left', 'skip': 'left'}},
     [1, -1, None, None, None, None], ['ref']),
])
def test_real_sources_absolute_arrows_exclusions_and_immutability(fixture, intent, arrows, flips):
    request, helper, context, resource_paths, response, paths = fixture
    bytes_before = [path.read_bytes() for path in paths]
    source_catalogs = [build_source_feature_catalog(SeqIO.read(path, 'genbank')) for path in paths]
    result = project(response, helper, intent)
    assert [row['afterArrow'] for row in result['records']] == arrows
    assert [row['recordKey'] for row in result['records'] if row['beforeReverseComplement'] != row['afterReverseComplement']] == flips
    assert [row['exclusion'] for row in result['records']] == [None, None, 'unknown_strand', 'skipped_by_user',
                                                             'skipped_no_candidate', 'skipped_unmappable']
    assert [path.read_bytes() for path in paths] == bytes_before
    assert [build_source_feature_catalog(SeqIO.read(path, 'genbank')) for path in paths] == source_catalogs
    assert set(response['plan']) == {'schema', 'groupId', 'reference', 'records'}
    # The final helper batch receives projected absolute directions, not a policy.
    final_context = {**context, 'orientations': {row['recordKey']: row['afterReverseComplement'] for row in result['records']}}
    final_response = resolve_similarity_alignment_payload(helper, projection=final_context, resource_paths=resource_paths)
    final_result = project(final_response, helper, intent, result['binding'])
    assert final_result['geometryValidated'] is True
    assert final_result['signature'] == result['signature']


def test_compound_crop_center_reference_fixed_all_y_and_idempotence(fixture):
    request, helper, context, resource_paths, response, _paths = fixture
    # Source envelope center is 180, despite the intron and retained crop fragment.
    ref_facts = response['projection']['records'][0]
    target_facts = response['projection']['records'][1]
    assert ref_facts['variants'][0]['anchors'][0]['displayCenter'] == 180
    assert ref_facts['variants'][1]['anchors'][0]['displayCenter'] == 820
    assert target_facts['variants'][1]['anchors'][0]['displayCenter'] == 170
    result = project(response, helper, {'mode': 'right'})
    assert result['reference']['beforeX'] == result['reference']['afterX']
    assert result['reference']['deltaX'] != 0
    candidate = _apply(request, result)
    candidate = replace(candidate, similarity_alignment=_plan_from_response(response))
    before_plan = plan_linear_request(request)
    after_plan = plan_linear_request(candidate)
    before_drawing = before_plan.build().drawing
    after_drawing = after_plan.build().drawing
    before_geometry = before_drawing._gbdraw_alignment_placements
    after_geometry = after_drawing._gbdraw_alignment_placements
    before_centers = project_similarity_alignment_centers(ResolvedRecordCollection(before_plan.records, before_plan.provenance), candidate.similarity_alignment)
    after_centers = after_plan.alignment_anchor_centers
    translations = _final_record_translations(record_keys=[record.record_key for record in candidate.records],
        placements=dict(enumerate(after_geometry)), layout=candidate.layout,
        similarity_alignment=candidate.similarity_alignment, anchor_centers=after_centers)
    # Keyed renderer metadata gives actual canvas geometry independently of the
    # projection; composition fitting is explicitly outside the logical frame.
    old_ref_x = (
        before_drawing._gbdraw_track_slot_geometry['records'][0]['axisXpx']
        - before_drawing._gbdraw_linear_composition_plan.placement_for('primary').dx
        + before_centers[0] * before_geometry[0].px_per_bp
    )
    assert result['reference']['beforeX'] == pytest.approx(old_ref_x)
    rendered_ref_x = (
        after_drawing._gbdraw_track_slot_geometry['records'][0]['axisXpx']
        - after_drawing._gbdraw_linear_composition_plan.placement_for('primary').dx
        + after_centers[0] * after_geometry[0].px_per_bp
    )
    assert rendered_ref_x == pytest.approx(old_ref_x)
    assert after_geometry[0].x_for_position(after_centers[0]) + translations[0][0] == pytest.approx(old_ref_x)
    assert after_geometry[1].x_for_position(after_centers[1]) + translations[1][0] == pytest.approx(old_ref_x)
    assert [placement.axis_y + xy[1] for placement, xy in zip(after_geometry, translations)] == pytest.approx(
        [placement.axis_y + translation.y for placement, translation in zip(before_geometry, request.layout.record_translations)])
    assert _final_record_translations(record_keys=[record.record_key for record in candidate.records],
        placements=dict(enumerate(after_geometry)), layout=candidate.layout,
        similarity_alignment=candidate.similarity_alignment, anchor_centers=after_centers) == translations
    # Ordinary Reverse continues to use current record orientation with the same plan.
    reversed_request = replace(candidate, records=(replace(candidate.records[0], presentation=replace(candidate.records[0].presentation, reverse_complement=True)),) + candidate.records[1:])
    reversed_plan = plan_linear_request(reversed_request)
    assert reversed_plan.alignment_anchor_centers[0] == 820
    assert reversed_request.similarity_alignment == candidate.similarity_alignment


def _plan_from_response(response):
    from gbdraw.layout.similarity_alignment import SimilarityAlignmentPlan, AlignmentRecordDecision, AlignmentDecisionStatus, AlignmentResolutionRationale
    raw = response['plan']
    return SimilarityAlignmentPlan(group_id=raw['groupId'], reference=_anchor(raw['reference'], 'ref'), records=tuple(
        AlignmentRecordDecision(row['recordKey'], AlignmentDecisionStatus(row['status']), AlignmentResolutionRationale(row['rationale']),
                                _anchor(row['anchor'], 'anchor') if row['anchor'] else None) for row in raw['records']))


def test_source_bytes_and_selector_changes_are_stale(fixture):
    _request, helper, context, resource_paths, response, paths = fixture
    old_binding = response['projection']['binding']
    record = SeqIO.read(paths[0], 'genbank')
    record.annotations['comment'] = 'Changed source bytes, same geometry.'
    SeqIO.write(record, paths[0], 'genbank')
    refreshed = resolve_similarity_alignment_payload(helper, projection=context, resource_paths=resource_paths)
    assert project(refreshed, helper, {'mode': 'keep'}, old_binding)['reason'] == 'binding_changed'
    changed = copy.deepcopy(context)
    changed['canonicalRequest']['records'][0]['selector'] = {'kind': 'recordId', 'value': 'ref'}
    selected = resolve_similarity_alignment_payload(helper, projection=changed, resource_paths=resource_paths)
    assert selected['projection']['binding'] != refreshed['projection']['binding']


def test_source_anchor_claims_and_canonical_coverage_rejected(fixture):
    _request, helper, context, resource_paths, _response, _paths = fixture
    helper = copy.deepcopy(helper)
    helper['members'][0]['sourceEnd'] = 261
    with pytest.raises(ValidationError, match='source anchor facts changed'):
        resolve_similarity_alignment_payload(helper, projection=context, resource_paths=resource_paths)
    helper['members'][0]['sourceEnd'] = 260
    helper['members'][0]['anchor']['stableFeatureSvgId'] = 'wrong'
    with pytest.raises(ValidationError, match='exactly one source feature'):
        resolve_similarity_alignment_payload(helper, projection=context, resource_paths=resource_paths)


def _refresh_fixture(request, helper, paths):
    helper = copy.deepcopy(helper)
    encoded = encode_canonical_request(request)
    for fact, record in zip(helper['records'], request.records, strict=True):
        fact['presentation']['reverseComplement'] = record.presentation.reverse_complement
        fact['region'] = None if record.region is None else dict(
            start=record.region.start, end=record.region.end, reverseComplement=record.region.reverse_complement)
    context = dict(canonicalRequest=encoded.payload, orientations=None)
    resources = {item.resource_id: str(item.source_path) for item in encoded.resources}
    response = resolve_similarity_alignment_payload(helper, projection=context, resource_paths=resources)
    return helper, context, resources, response


def test_real_minority_reference_and_already_reversed_absolute_direction(fixture):
    request, helper, _context, _resources, _response, paths = fixture
    target = replace(request.records[1], region=replace(request.records[1].region, reverse_complement=False))
    request = replace(request, records=(request.records[0], target) + request.records[2:])
    helper, _context, _resources, response = _refresh_fixture(request, helper, paths)
    right = project(response, helper, {'mode': 'right'})
    assert [row['beforeArrow'] for row in right['records'][:2]] == [-1, 1]
    assert [row['recordKey'] for row in right['records'] if row['beforeReverseComplement'] != row['afterReverseComplement']] == ['ref']
    already = _apply(request, right)
    helper, _context, _resources, response = _refresh_fixture(already, helper, paths)
    unchanged = project(response, helper, {'mode': 'right'})
    assert unchanged['reference']['deltaX'] == 0
    assert all(row['beforeReverseComplement'] == row['afterReverseComplement'] for row in unchanged['records'])
    left_ref = project(response, helper, {'mode': 'custom', 'byRecordKey': {'ref': 'left'}})
    assert left_ref['records'][0]['beforeReverseComplement'] is True
    assert left_ref['records'][0]['afterReverseComplement'] is False
    assert [row['afterArrow'] for row in left_ref['records'][:2]] == [-1, 1]


def test_cropped_compound_identity_uses_full_source_center_and_cosmetics_do_not_stale(fixture):
    request, helper, _context, _resources, _response, paths = fixture
    reference = replace(request.records[0], region=parse_region_spec('ref:151-350'))
    request = replace(request, records=(reference,) + request.records[1:])
    helper, context, resources, response = _refresh_fixture(request, helper, paths)
    record = response['projection']['records'][0]
    assert record['variants'][0]['anchors'][0]['displayCenter'] == 30
    assert record['variants'][1]['anchors'][0]['displayCenter'] == 170
    result = project(response, helper, {'mode': 'right'})
    assert result['reference']['beforeX'] == result['reference']['afterX']
    assert result['records'][0]['anchor'] == helper['reference']
    cosmetic = copy.deepcopy(context)
    cosmetic['canonicalRequest']['records'][0]['presentation']['label'] = 'Changed display label'
    cosmetic['canonicalRequest']['records'][0]['presentation']['subtitle'] = 'Cosmetic subtitle'
    refreshed = resolve_similarity_alignment_payload(helper, projection=cosmetic, resource_paths=resources)
    assert refreshed['projection']['binding'] == response['projection']['binding']
    assert project(refreshed, helper, {'mode': 'keep'}, result['binding'])['status'] == 'projected'


def test_realign_materialized_geometry_has_no_accumulated_offset(fixture):
    request, helper, _context, _resources, response, paths = fixture
    intent = {'mode': 'right'}
    preview = project(response, helper, intent)
    context = dict(canonicalRequest=encode_canonical_request(request).payload,
                   orientations={row['recordKey']: row['afterReverseComplement'] for row in preview['records']})
    resources = {item.resource_id: str(item.source_path) for item in encode_canonical_request(request).resources}
    final_response = resolve_similarity_alignment_payload(helper, projection=context, resource_paths=resources)
    result = project(final_response, helper, intent)
    candidate = replace(_apply(request, result), similarity_alignment=_plan_from_response(final_response))
    plan = plan_linear_request(candidate)
    geometry = plan.build().drawing._gbdraw_alignment_placements
    absolute = _final_record_translations(record_keys=[record.record_key for record in candidate.records],
        placements=dict(enumerate(geometry)), layout=candidate.layout,
        similarity_alignment=candidate.similarity_alignment, anchor_centers=plan.alignment_anchor_centers)
    # The existing bridge materializes these rendered translations on the next Align.
    materialized = replace(candidate, similarity_alignment=None,
        layout=replace(candidate.layout, record_translations=tuple(
            LinearRecordTranslation(record.record_key, *xy) for record, xy in zip(candidate.records, absolute, strict=True))))
    helper, _context, _resources, response = _refresh_fixture(materialized, helper, paths)
    repeated = project(response, helper, {'mode': 'keep'})
    assert repeated['reference']['deltaX'] == 0
    second = replace(_apply(materialized, repeated), similarity_alignment=_plan_from_response(response))
    second_plan = plan_linear_request(second)
    second_geometry = second_plan.build().drawing._gbdraw_alignment_placements
    actual = _final_record_translations(record_keys=[record.record_key for record in second.records],
        placements=dict(enumerate(second_geometry)), layout=second.layout,
        similarity_alignment=second.similarity_alignment, anchor_centers=second_plan.alignment_anchor_centers)
    assert actual == absolute


def test_concurrent_baseline_placement_change_is_explicitly_stale(fixture):
    _request, helper, context, resources, response, _paths = fixture
    changed = copy.deepcopy(context)
    changed['canonicalRequest']['layout']['recordTranslations'][0]['x'] += 17
    refreshed = resolve_similarity_alignment_payload(helper, projection=changed, resource_paths=resources)
    stale = project(refreshed, helper, {'mode': 'right'}, response['projection']['binding'])
    assert stale['status'] == 'stale' and stale['reason'] == 'binding_changed'
