"""Evaluate browser rule drafts through the renderer's Python rule owners."""
from __future__ import annotations

import json
from types import SimpleNamespace

from pandas import DataFrame

from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.selector_values import iter_specific_color_rules
from gbdraw.labels.filtering import _build_label_override_rules, _resolve_label_override


def evaluate_rules_json(features_json: str, rules_json: str, kind: str = "color") -> str:
    features = json.loads(features_json)
    rules = json.loads(rules_json)
    if kind == "color":
        rows = [dict(feature_type=r["feat"], qualifier_key=r["qual"], value=r["val"],
                     color=str(i), caption="") for i, r in enumerate(rules)]
        defaults = DataFrame(columns=["feature_type", "color"])
        compiled, _ = preprocess_color_tables(DataFrame(rows), defaults)
    elif kind == "label":
        rows = [dict(record_id=r["recordId"], feature_type=r["featureType"],
                     qualifier=r["qualifier"], value=r["valueRegex"], label_text=str(i))
                for i, r in enumerate(rules)]
        compiled = _build_label_override_rules(DataFrame(rows)) or []
    else:
        raise ValueError("Unknown rule kind")
    matches, winners, priorities = [], [], []
    for item in features:
        feature = SimpleNamespace(type=item["type"], qualifiers=item["qualifiers"])
        selector = item["selector"]
        record = item.get("record", "")
        if kind == "color":
            ordered = list(iter_specific_color_rules(feature, compiled, record, selector=selector))
            matched = list(dict.fromkeys(int(match[0]) for match in ordered))
            ranks = [None] * len(rules)
            for color, _, rank in ordered:
                if ranks[int(color)] is None:
                    ranks[int(color)] = rank
            priorities.append(ranks)
            matches.append(matched)
            winners.append(matched[0] if matched else -1)
        else:
            winner = _resolve_label_override(feature, feature.type, feature.qualifiers,
                                             item.get("label", ""), compiled, record, selector=selector)
            winners.append(int(winner) if winner is not None else -1)
    return json.dumps({"matches": matches, "winners": winners, "priorities": priorities})
