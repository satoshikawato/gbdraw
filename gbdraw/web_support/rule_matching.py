"""Evaluate browser rule drafts through the renderer's Python rule owners."""
from __future__ import annotations

import json
from types import SimpleNamespace

from pandas import DataFrame

from gbdraw.features.colors import normalize_specific_color_captions, preprocess_color_tables
from gbdraw.features.selector_values import iter_specific_color_rules
from gbdraw.features.visibility import (
    _first_matching_visibility_rule,
    compile_feature_visibility_rules,
)
from gbdraw.labels.filtering import _build_label_override_rules, _resolve_label_override


def _visibility_table(rules: list[dict]) -> DataFrame:
    return DataFrame([dict(
        record_id=rule["recordId"], feature_type=rule["featureType"],
        qualifier=rule["qualifier"], value=rule["value"], action=rule["action"],
    ) for rule in rules])


def _compile_visibility_rows(rules: list[dict]) -> list[dict | None]:
    """The draft rows as Generate compiles the visibility table, with None for
    a row that Generate skips as a header. A table Generate rejects raises
    Generate's error, which names the row."""
    compile_feature_visibility_rules(_visibility_table(rules))
    rows = (compile_feature_visibility_rules(_visibility_table([rule])) for rule in rules)
    return [compiled[0] if compiled else None for compiled in rows]


def evaluate_rules_json(features_json: str, rules_json: str, kind: str = "color") -> str:
    features = json.loads(features_json)
    rules = json.loads(rules_json)
    if kind == "color-captions":
        table = DataFrame(
            [
                dict(
                    feature_type=r["feat"],
                    qualifier_key=r["qual"],
                    value=r["val"],
                    color=r["color"],
                    caption=r.get("cap", ""),
                )
                for r in rules
            ]
        )
        normalized = normalize_specific_color_captions(table)
        assert normalized is not None  # a DataFrame input is never None
        captions = list(normalized["caption"]) if len(rules) else []
        return json.dumps(
            {
                "rules": [
                    dict(rule, cap=caption) for rule, caption in zip(rules, captions)
                ]
            }
        )
    if kind == "color":
        rows = [
            dict(
                feature_type=r["feat"],
                qualifier_key=r["qual"],
                value=r["val"],
                color=str(i),
                caption="",
            )
            for i, r in enumerate(rules)
        ]
        defaults = DataFrame(columns=["feature_type", "color"])
        color_rules, _ = preprocess_color_tables(DataFrame(rows), defaults)
    elif kind == "visibility":
        visibility_rules = _compile_visibility_rows(rules)
    elif kind == "label":
        rows = [dict(record_id=r["recordId"], feature_type=r["featureType"],
                     qualifier=r["qualifier"], value=r["valueRegex"], label_text=str(i))
                for i, r in enumerate(rules)]
        label_rules = _build_label_override_rules(DataFrame(rows)) or []
    else:
        raise ValueError("Unknown rule kind")
    matches, winners, priorities = [], [], []
    for item in features:
        feature = SimpleNamespace(type=item["type"], qualifiers=item["qualifiers"])
        selector = item["selector"]
        record = item.get("record", "")
        if kind == "color":
            ordered = list(iter_specific_color_rules(feature, color_rules, record, selector=selector))
            matched = list(dict.fromkeys(int(match[0]) for match in ordered))
            ranks = [None] * len(rules)
            for color, _, rank in ordered:
                if ranks[int(color)] is None:
                    ranks[int(color)] = rank
            priorities.append(ranks)
            matches.append(matched)
            winners.append(matched[0] if matched else -1)
        elif kind == "visibility":
            # Every rule this feature matches; the first in table order decides.
            matches.append([
                index for index, rule in enumerate(visibility_rules)
                if rule is not None
                and _first_matching_visibility_rule(feature, [rule], record, selector=selector)
            ])
        else:
            winner = _resolve_label_override(feature, feature.type, feature.qualifiers,
                                             item.get("label", ""), label_rules, record, selector=selector)
            winners.append(int(winner) if winner is not None else -1)
    return json.dumps({"matches": matches, "winners": winners, "priorities": priorities})
