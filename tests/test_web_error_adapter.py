"""Native failure identity, bounded correction facts, and browser privacy."""
from __future__ import annotations

import json
import logging
from pathlib import Path
import re
from types import SimpleNamespace

from pandas import DataFrame
import pytest

from gbdraw.api import DepthTrackInput
from gbdraw.exceptions import ComparisonIdentityError, ParseError, ValidationError
from gbdraw.features.visibility import compile_feature_visibility_rules
from gbdraw.io.record_select import parse_record_selector, select_record
from gbdraw.labels.filtering import _build_label_override_rules
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup
from gbdraw.web_support.error_adapter import serialize_web_error, private_web_execution


@pytest.fixture(scope='module')
def helpers():
    source = (Path(__file__).resolve().parents[1] / 'gbdraw/web/js/app/python-helpers.js').read_text()
    namespace = {}
    exec(source.split('`', 1)[1].rsplit('`', 1)[0], namespace)
    return namespace


@pytest.mark.parametrize(('value', 'indexes', 'reason'), [
    ('private;', '', 'EMPTY_ENDPOINT'), ('private;second', '0', 'INDEX_ALIGNMENT'),
    ('private', '-1', 'SOURCE_INDEX'), ('private', 'invalid', 'SOURCE_INDEX'),
    ('private', '0', 'SOURCE_VIEW_CONFLICT'),
])
def test_comparison_identity_producer_preserves_valueerror(value, indexes, reason):
    group = object.__new__(PairWiseMatchGroup)
    group.records = [object()]
    group.feature_dom_index = SimpleNamespace(by_source_index={(0, 0): 'other'}, by_view_id={(0, 'private'): ['private']})
    with pytest.raises(ValueError) as caught:
        group._rendered_feature_svg_id_values(value, record_index=0, feature_index_value=indexes)
    assert isinstance(caught.value, ComparisonIdentityError)
    payload = serialize_web_error(caught.value, operation='generate', stage='render')
    assert payload['code'] == 'COMPARISON_IDENTITY'
    assert payload['context'] == {'reason': reason}
    assert 'private' not in json.dumps(payload)


@pytest.mark.parametrize('catalog', [[], [{'type': 'CDS', 'qualifiers': {'product': ['unrelated']}, 'selector': {}, 'record': 'private'}]])
@pytest.mark.parametrize('kind', ['color', 'label'])
def test_actual_helper_retains_explicit_regex_cause(helpers, kind, catalog, capsys):
    pattern = '😀[PRIVATE_PATTERN_SENTINEL'
    rule = dict(feat='CDS', qual='product', val=pattern) if kind == 'color' else dict(
        recordId='*', featureType='CDS', qualifier='product', valueRegex=pattern)
    # The native function still raises its original re.error/ParseError with cause.
    with pytest.raises((re.error, ParseError)) as native:
        helpers['evaluate_rules_json'](json.dumps(catalog), json.dumps([rule]), kind)
    if kind == 'label':
        assert isinstance(native.value.__cause__, re.error)
    capsys.readouterr()
    bounded = json.loads(helpers['call_web_json_helper']('evaluate_rules_json', json.dumps(catalog), json.dumps([rule]), kind))['error']
    assert bounded['code'] == 'REGEX_SYNTAX'
    assert bounded['operation'] == 'evaluateRules'
    assert bounded['stage'] == 'rule-validation'
    assert bounded['context']['position'] == 1  # Python characters, not UTF-16
    assert bounded['context']['reason'] == 'UNTERMINATED_SET'
    if kind == 'label':
        assert bounded['context']['row'] == 1
    assert 'PRIVATE_PATTERN_SENTINEL' not in json.dumps(bounded)
    assert capsys.readouterr() == ('', '')


@pytest.mark.parametrize(('call', 'code', 'context'), [
    (lambda: parse_record_selector('#0'), 'RECORD_SELECTION', {'reason': 'POSITIVE_INTEGER'}),
    (lambda: select_record([], parse_record_selector('#2')), 'RECORD_SELECTION', {'reason': 'OUT_OF_RANGE', 'recordCount': 0}),
    (lambda: select_record([SimpleNamespace(id='private'), SimpleNamespace(id='private')], parse_record_selector('private')), 'RECORD_SELECTION', {'reason': 'AMBIGUOUS'}),
    (lambda: DepthTrackInput(source=[]), 'DEPTH_INVALID', {'field': 'source', 'reason': 'REQUIRED'}),
    (lambda: DepthTrackInput(source='file', height=-1), 'INPUT_INVALID', {'field': 'height', 'reason': 'POSITIVE_OR_AUTO'}),
    (lambda: _build_label_override_rules(DataFrame([dict(record_id='*', feature_type='CDS', qualifier='', value='valid', label_text='private')])), 'TABLE_INVALID', {'field': 'qualifier', 'row': 1, 'reason': 'REQUIRED'}),
])
def test_known_native_validation_correction_is_not_unknown(call, code, context):
    with pytest.raises((ValueError, ValidationError)) as caught:
        call()
    result = serialize_web_error(caught.value, operation='generate', stage='request-validation')
    assert result['code'] == code
    assert result['context'] == context


@pytest.mark.parametrize(('message', 'code', 'context'), [
    ('min_depth must be <= max_depth.', 'INPUT_INVALID', {'field': 'min_depth', 'reason': 'ORDER'}),
    ('bitscore must be a finite value >= 0', 'INPUT_INVALID', {'field': 'bitscore', 'reason': 'NONNEGATIVE'}),
    ('Malformed line in label override file \'PRIVATE_PATH_SENTINEL\' at line 4: expected 5 columns.', 'TABLE_INVALID', {'row': 4, 'columnCount': 5}),
])
def test_known_template_adapter_keeps_only_allowlisted_facts(message, code, context):
    result = serialize_web_error(ValidationError(message), operation='generate', stage='render')
    assert result['code'] == code
    assert result['context'] == context
    assert 'PRIVATE_' not in json.dumps(result)


def test_unknown_and_unmapped_validation_are_distinct_and_private(helpers, capsys):
    sentinel = 'PRIVATE_SENTINEL' * 5000
    error = RuntimeError('unknown extension ' + sentinel)
    error.__cause__ = RuntimeError(sentinel)
    error._web_error_secondary = [{'code': 'CLEANUP_FAILED', 'stage': 'cleanup', 'message': sentinel}]
    assert serialize_web_error(error, operation='generate', stage='render') == {
        'code': 'UNKNOWN', 'operation': 'generate', 'stage': 'render', 'context': {}, 'secondary': []}
    assert serialize_web_error(ValidationError(sentinel), operation='generate', stage='render')['code'] == 'VALIDATION_UNCLASSIFIED'
    def fail():
        print(sentinel)
        logging.error(sentinel)
        raise error
    helpers['_WEB_JSON_HELPERS']['measure_legend_text_json'] = (fail, 'measureLegendText')
    result = json.loads(helpers['call_web_json_helper']('measure_legend_text_json'))
    assert result['error']['code'] == 'UNKNOWN'
    assert sentinel not in json.dumps(result)
    assert capsys.readouterr() == ('', '')
    # The scope restores native logging and does not change native cause/capture.
    previous = logging.root.manager.disable
    with pytest.raises(RuntimeError), private_web_execution():
        raise error
    assert logging.root.manager.disable == previous
    assert isinstance(error.__cause__, RuntimeError)


def test_visibility_syntax_cause_and_missing_position_remain_native():
    with pytest.raises(ParseError) as caught:
        compile_feature_visibility_rules(DataFrame([dict(record_id='*', feature_type='CDS', qualifier='product', value='[', action='show')]))
    assert isinstance(caught.value.__cause__, re.error)
    result = serialize_web_error(caught.value, operation='generate', stage='render')
    assert result['code'] == 'REGEX_SYNTAX'
    assert result['context']['row'] == 1
    unknown_position = serialize_web_error(re.error('synthetic private reason'), operation='evaluateRules', stage='helper')
    assert 'position' not in unknown_position['context']
    assert unknown_position['context']['reason'] == 'SYNTAX_ERROR'


def test_native_visibility_action_retains_row_and_correction():
    with pytest.raises(ValidationError) as caught:
        compile_feature_visibility_rules(DataFrame([dict(record_id='*', feature_type='CDS', qualifier='product', value='valid', action='PRIVATE_ACTION')]))
    model = serialize_web_error(caught.value, operation='generate', stage='render')
    assert model['code'] == 'TABLE_INVALID'
    assert model['context'] == {'field': 'action', 'row': 1, 'reason': 'VISIBILITY_ACTION'}
    assert 'PRIVATE_' not in json.dumps(model)


def test_render_json_failure_reports_actual_request_stage(helpers, capsys):
    model = helpers['run_canonical_request_wrapper']('{PRIVATE_JSON', '{}', '/PRIVATE_PATH')['error']
    assert model == {'code': 'HELPER_PROTOCOL', 'operation': 'generate',
                     'stage': 'request-validation', 'context': {'reason': 'JSON_FORMAT'}}
    assert capsys.readouterr() == ('', '')


def test_explicit_cause_cycles_and_invalid_identifiers_stay_bounded():
    error = RuntimeError('PRIVATE_ERROR')
    error.__cause__ = error
    error._web_error_stage = ['PRIVATE_STAGE']
    model = serialize_web_error(error, operation=['PRIVATE_OPERATION'], stage=['PRIVATE_STAGE'])
    assert model == {'code': 'UNKNOWN', 'operation': 'unknown', 'stage': 'unknown', 'context': {}}


@pytest.mark.parametrize(('validator', 'field', 'reason'), [
    ('normalize_collinearity_anchor_mode', 'collinear_anchor_mode', 'COLLINEAR_ANCHOR_MODE'),
    ('normalize_collinearity_color_mode', 'collinear_color_mode', 'COLLINEAR_COLOR_MODE'),
    ('normalize_collinearity_search_scope', 'collinear_search_scope', 'ADJACENT_ALL'),
])
def test_native_comparison_choice_corrections(validator, field, reason):
    from gbdraw.analysis import collinearity
    with pytest.raises(ValidationError) as caught:
        getattr(collinearity, validator)('PRIVATE_OPTION')
    result = serialize_web_error(caught.value, operation='convertLosatpPairsToGenomicPayload', stage='helper')
    assert result['code'] == 'INPUT_INVALID'
    assert result['context'] == {'field': field, 'reason': reason}


@pytest.mark.parametrize(('region', 'reason'), [
    ({'start': 5, 'end': 1, 'reverseComplement': False}, 'ORDER'),
    ({'start': 1, 'end': 11, 'reverseComplement': False}, 'RECORD_BOUNDS'),
])
def test_native_alignment_region_correction(region, reason):
    from gbdraw.web_support.similarity_alignment import _record_fact
    with pytest.raises(ValidationError) as caught:
        _record_fact({'recordKey': 'PRIVATE_RECORD', 'recordLength': 10, 'region': region,
                      'presentation': {'reverseComplement': False}}, 'records[0]')
    result = serialize_web_error(caught.value, operation='resolveSimilarityAlignment', stage='helper')
    assert result['code'] == 'INPUT_INVALID'
    assert result['context'] == {'field': 'region', 'reason': reason}


def test_missing_gff_fasta_match_keeps_correction_without_record_name():
    error = ParseError(
        "No matching FASTA record found for GFF record PRIVATE_RECORD. "
        "Please ensure that all GFF records have corresponding FASTA entries."
    )
    model = serialize_web_error(error, operation='generate', stage='render')
    assert model['code'] == 'FASTA_REQUIRED'
    assert model['context'] == {'reason': 'GFF_FASTA_MATCH'}
    assert 'PRIVATE_' not in json.dumps(model)


def test_unknown_config_path_keeps_correction_without_private_path():
    from gbdraw.config.modify import validate_config_overrides
    with pytest.raises(ValidationError) as caught:
        validate_config_overrides({'PRIVATE_PATH': 1})
    model = serialize_web_error(caught.value, operation='validateConfigOverrides', stage='helper')
    assert model['code'] == 'INPUT_INVALID'
    assert model['context'] == {'field': 'configOverrides', 'reason': 'UNKNOWN_CONFIG_PATH'}
    assert 'PRIVATE_' not in json.dumps(model)
