"""Read-only synthetic calls to existing Python regex owners; no render or UI."""
import contextlib
import io
import json
import logging
import re

from Bio.SeqFeature import SeqFeature, SimpleLocation
from pandas import DataFrame

from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.selector_values import build_feature_selector_values, find_specific_color_rule
from gbdraw.features.visibility import compile_feature_visibility_rules, _first_matching_visibility_rule
from gbdraw.labels.filtering import _build_whitelist_map, _matches_any_pattern, _build_label_override_rules, _resolve_label_override
from gbdraw.web_support.rule_matching import evaluate_rules_json

PATTERNS = ["(?i)NADH", "(?P<enzyme>NADH)", "[", "(?<enzyme>NADH)", r"NADH\Z", r"\bβ", "i", "β"]
CATALOGS = {"corpus": ["NADH", "nadh", "β-lactamase", "ı", "i", "İ", "unrelated"],
            "unrelated": ["unrelated"], "empty": []}


def capture(call):
    # Capture only synthetic producer logging, so evidence records exposure without printing it.
    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    root = logging.getLogger()
    previous = root.handlers[:]
    root.handlers = [handler]
    try:
        with contextlib.redirect_stderr(stream):
            result = call()
        return {"accepted": True, "targets": result}
    except Exception as error:
        cause = error if isinstance(error, re.error) else error.__cause__
        return {"accepted": False, "type": type(error).__name__,
                "causeType": type(cause).__name__ if cause is not None else None,
                "pythonCharacterPosition": getattr(cause, "pos", None),
                "producerLoggedSyntheticPattern": any(p in stream.getvalue() for p in ["[", "(?<enzyme>NADH)"])}
    finally:
        root.handlers = previous


rows = []
for pattern in PATTERNS:
    for catalog, texts in CATALOGS.items():
        features = [SeqFeature(SimpleLocation(i * 100, i * 100 + 50, strand=1), type="CDS",
                               qualifiers={"product": [text]}) for i, text in enumerate(texts)]
        payload = [dict(type=f.type, qualifiers=f.qualifiers, record="synthetic", label=text,
                        selector=build_feature_selector_values(f, "synthetic")) for f, text in zip(features, texts)]
        color = [{"feat": "CDS", "qual": "product", "val": pattern}]
        label = [{"recordId": "*", "featureType": "CDS", "qualifier": "product", "valueRegex": pattern}]

        def native_color():
            compiled, _ = preprocess_color_tables(DataFrame([dict(feature_type="CDS", qualifier_key="product", value=pattern, color="hit")]),
                                                   DataFrame(columns=["feature_type", "color"]))
            return [i for i, f in enumerate(features) if find_specific_color_rule(f, compiled, "synthetic")]

        def native_label():
            compiled = _build_label_override_rules(DataFrame([dict(record_id="*", feature_type="CDS", qualifier="product", value=pattern, label_text="hit")]))
            return [i for i, (f, text) in enumerate(zip(features, texts)) if _resolve_label_override(f, f.type, f.qualifiers, text, compiled, "synthetic") is not None]

        def whitelist():
            compiled = _build_whitelist_map(DataFrame([dict(feature_type="CDS", qualifier="product", keyword=pattern)]))
            return [i for i, text in enumerate(texts) if _matches_any_pattern(text, compiled["CDS"]["product"])]

        def visibility():
            compiled = compile_feature_visibility_rules(DataFrame([dict(record_id="*", feature_type="CDS", qualifier="product", value=pattern, action="hide")]))
            return [i for i, f in enumerate(features) if _first_matching_visibility_rule(f, compiled, "synthetic")]

        result = {"pattern": pattern, "catalog": catalog}
        for kind, rules in [("color", color), ("label", label)]:
            result[kind] = capture(lambda: [i for i, winner in enumerate(json.loads(evaluate_rules_json(json.dumps(payload), json.dumps(rules), kind))["winners"]) if winner >= 0])
        result["nativeColor"] = capture(native_color)
        result["nativeLabel"] = capture(native_label)
        result["whitelist"] = capture(whitelist)
        result["visibility"] = capture(visibility)
        for kind, native in [("color", "nativeColor"), ("label", "nativeLabel")]:
            assert result[kind]["accepted"] == result[native]["accepted"]
            if result[kind]["accepted"]:
                assert result[kind]["targets"] == result[native]["targets"]
        assert all(result[k]["accepted"] == (pattern not in ["[", "(?<enzyme>NADH)"])
                   for k in ["color", "label", "nativeColor", "nativeLabel", "whitelist", "visibility"])
        rows.append(result)

print(json.dumps({"scope": "native/helper regex owners only; no Worker, browser, Generate, or original-audit reproduction", "corpus": CATALOGS, "rows": rows}, ensure_ascii=False, indent=2))
