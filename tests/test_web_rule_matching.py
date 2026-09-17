"""Browser drafts and native rendering must use the same regex engine and selectors."""
import json
import re

import pytest
from Bio.SeqFeature import SeqFeature, SimpleLocation
from pandas import DataFrame

from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.selector_values import build_feature_selector_values, find_specific_color_rule
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
