"""Browser drafts and native rendering must use the same regex engine and selectors."""
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.SeqFeature import SeqFeature, SimpleLocation
from pandas import DataFrame

from gbdraw.exceptions import ParseError
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.selector_values import build_feature_selector_values, find_specific_color_rule
from gbdraw.features.visibility import (
    _first_matching_visibility_rule,
    compile_feature_visibility_rules,
    should_render_feature,
)
from gbdraw.labels.filtering import _build_label_override_rules, _resolve_label_override
from gbdraw.web_support.rule_matching import evaluate_rules_json

CORPUS = ['NADH', 'NADH dehydrogenase', 'β-lactamase', 'ı', 'i', 'İ', 'unrelated']
PATTERNS = ['(?i)NADH', '(?P<enzyme>NADH)', r'NADH\Z', r'\bβ', 'i', r'^NADH$', '[']


@pytest.mark.parametrize('pattern', PATTERNS + ['(?<enzyme>NADH)'])
def test_browser_rule_targets_equal_native_color_and_label_rules(pattern):
    features = [SeqFeature(SimpleLocation(i * 100, i * 100 + 50, strand=1), type='CDS',
                           qualifiers={'product': [text]}) for i, text in enumerate(CORPUS)]
    payload = [dict(type=f.type, qualifiers=f.qualifiers,
                    selector=build_feature_selector_values(f, 'rec'), record='rec', label=text)
               for f, text in zip(features, CORPUS)]
    color = [dict(feat='CDS', qual='product', val=pattern)]
    label = [dict(recordId='*', featureType='CDS', qualifier='product', valueRegex=pattern)]
    try:
        re.compile(pattern, re.I)
    except re.error:
        for kind, rules in [('color', color), ('label', label)]:
            with pytest.raises(Exception, match='(unterminated|unknown extension)'):
                evaluate_rules_json(json.dumps(payload), json.dumps(rules), kind)
        return
    color_map, _ = preprocess_color_tables(
        DataFrame([dict(feature_type='CDS', qualifier_key='product', value=pattern, color='hit')]),
        DataFrame([dict(feature_type='CDS', color='miss')]))
    label_rules = _build_label_override_rules(DataFrame([
        dict(record_id='*', feature_type='CDS', qualifier='product', value=pattern, label_text='hit')]))
    color_expected = [0 if find_specific_color_rule(f, color_map, 'rec') else -1 for f in features]
    label_expected = [0 if _resolve_label_override(f, f.type, f.qualifiers, text, label_rules, 'rec') else -1
                      for f, text in zip(features, CORPUS)]
    assert json.loads(evaluate_rules_json(json.dumps(payload), json.dumps(color)))['winners'] == color_expected
    assert json.loads(evaluate_rules_json(json.dumps(payload), json.dumps(label), 'label'))['winners'] == label_expected
    assert color_expected == label_expected


def test_empty_catalog_still_validates_all_rules():
    with pytest.raises(re.error):
        evaluate_rules_json('[]', '[{"feat":"CDS","qual":"product","val":"(?<enzyme>NADH)"}]')


# R4: tests/web/feature-drawn-resolver.test.mjs runs the same cases through the
# JavaScript resolver of the Label On dialog and the live preview.
_DRAWN = json.loads(
    (Path(__file__).parent / "fixtures" / "feature_drawn_cases.json").read_text(encoding="utf-8")
)


def _drawn_feature(name):
    spec = _DRAWN["features"][name]
    start, end, strand = spec["location"]
    feature = SeqFeature(
        SimpleLocation(start, end, strand=strand), type=spec["type"], qualifiers=spec["qualifiers"]
    )
    return feature, spec


def _visibility_rules(case):
    return compile_feature_visibility_rules(DataFrame(
        [dict(record_id=r["recordId"], feature_type=r["featureType"], qualifier=r["qualifier"],
              value=r["value"], action=r["action"]) for r in case.get("visibilityRules", [])],
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    ))


@pytest.mark.parametrize("name", sorted(_DRAWN["features"]))
def test_shared_drawn_vectors_hold_the_drawn_selector_values(name):
    feature, spec = _drawn_feature(name)
    values = build_feature_selector_values(feature, spec["recordId"])
    assert spec["drawnSelector"] == {
        "hash": values.get("hash"),
        "location": values.get("location"),
        "recordLocation": values.get("record_location"),
    }


@pytest.mark.parametrize("case", _DRAWN["cases"], ids=lambda case: case["name"])
def test_should_render_feature_answers_the_shared_drawn_vectors(case):
    feature, spec = _drawn_feature(case["feature"])
    color_rules = case.get("colorRules", [])
    color_map = preprocess_color_tables(
        DataFrame([dict(feature_type=r["feat"], qualifier_key=r["qual"], value=r["val"],
                        color="#123456", caption="") for r in color_rules]),
        DataFrame(columns=["feature_type", "color"]),
    )[0] if color_rules else None
    override = SimpleNamespace(feature_visibility=case["override"]) if case.get("override") else None
    assert should_render_feature(
        feature,
        case["selectedFeatures"],
        _visibility_rules(case),
        record_id=spec["recordId"],
        specific_color_rules=color_map,
        feature_override=override,
    ) is case["drawn"]


@pytest.mark.parametrize(
    "case", [case for case in _DRAWN["cases"] if case.get("visibilityRules")], ids=lambda case: case["name"]
)
def test_visibility_helper_matches_the_rules_generate_matches(case):
    """The helper reads the drawn selector values the catalog gives the Web."""
    feature, spec = _drawn_feature(case["feature"])
    drawn = spec["drawnSelector"]
    payload = [dict(type=spec["type"], qualifiers=spec["qualifiers"], record=spec["recordId"], selector={
        "hash": drawn["hash"], "location": drawn["location"], "record_location": drawn["recordLocation"],
    })]
    rules = case["visibilityRules"]
    matches = json.loads(evaluate_rules_json(json.dumps(payload), json.dumps(rules), "visibility"))["matches"][0]
    compiled = _visibility_rules(case)
    native = [index for index, rule in enumerate(compiled)
              if _first_matching_visibility_rule(feature, [rule], spec["recordId"])]
    assert matches == native


def test_visibility_helper_rejects_a_table_with_the_error_and_row_generate_reports():
    """OV-19: a rule edit reports what Generate reports for the same table."""
    rules = [
        dict(recordId="*", featureType="CDS", qualifier="locus_tag", value="^fl", action="off"),
        dict(recordId="*", featureType="CDS", qualifier="locus_tag", value="^fl1(", action="show"),
    ]
    table = DataFrame([dict(record_id=r["recordId"], feature_type=r["featureType"], qualifier=r["qualifier"],
                            value=r["value"], action=r["action"]) for r in rules])
    with pytest.raises(ParseError, match="at row 2:") as generate:
        compile_feature_visibility_rules(table)
    with pytest.raises(ParseError) as helper:
        evaluate_rules_json("[]", json.dumps(rules), "visibility")
    assert str(helper.value) == str(generate.value)
