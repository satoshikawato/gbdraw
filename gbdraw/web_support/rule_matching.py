"""Evaluate browser rule drafts through the renderer's Python rule owners."""
from __future__ import annotations

import json
import re
from typing import Any, Callable, Literal, Sequence, TypedDict, is_typeddict

from types import SimpleNamespace

from pandas import DataFrame

from gbdraw.contract_shapes import Shape, ref, typed_dict_shape, violations
from gbdraw.features.colors import build_specific_color_rules
from gbdraw.features.selector_values import iter_specific_color_rules
from gbdraw.features.visibility import (
    _first_matching_visibility_rule,
    compile_feature_visibility_rules,
)
from gbdraw.labels.filtering import DEFAULT_LABEL_PRIORITY, _build_label_override_rules, _resolve_label_override


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


class RuleFeatureSelector(TypedDict):
    """The drawn values a `hash`, `location`, or `record_location` rule reads."""

    hash: str | None
    location: str | None
    record_location: str | None


class RuleFeaturePayload(TypedDict):
    """A feature as Python's rule evaluation reads it (gbdraw/web/js/app/rule-matching.js `ruleFeaturePayload`)."""

    type: str
    record: str
    qualifiers: dict[str, list[str]]
    selector: RuleFeatureSelector
    label: str


class ColorRuleCondition(TypedDict):
    """What a Specific color rule matches: a feature type, a qualifier, and a pattern."""

    feat: str
    qual: str
    val: str


class ColorRuleRow(ColorRuleCondition):
    """A Specific color rule row as the Web's rule table holds it."""

    color: str
    cap: str


class SynthesizeColorRulesRequest(TypedDict):
    """One color edit: color the `selected` features of the catalog `features`.

    `scope` "group" colors them as one group (exact conditions where they fit,
    Owner R-07); "features" colors each by itself (a `hash` rule per feature, as
    "This feature only" does). `labels` holds, for each selected feature in the
    same order, the label texts it shows (its displayed and its annotation
    label). `rules` is the rule table at the write's turn, in order.
    """

    scope: Literal["group", "features"]
    features: list[RuleFeaturePayload]
    selected: list[int]
    labels: list[list[str]]
    rules: list[ColorRuleCondition]
    color: str
    caption: str


class SynthesizedColorRules(TypedDict):
    """The next rule table: a number keeps that row of the request's `rules`, an
    object is a new row. `declined` lists the selected features no rule can
    color (no exact condition and no hash)."""

    table: list[int | ColorRuleRow]
    exactRules: int
    hashRules: int
    declined: list[int]


def describe_rule_synthesis() -> dict[str, Any]:
    """The request and the answer of :func:`synthesize_color_rules_json` as named
    shapes of :mod:`gbdraw.contract_shapes` (tools/generate_web_contracts.py)."""

    types_: dict[str, Shape] = {}

    def special(hint: object) -> Shape | None:
        if isinstance(hint, type) and is_typeddict(hint):
            if hint.__name__ not in types_:
                types_[hint.__name__] = {"kind": "unknown"}  # reserves its place before its members'
                types_[hint.__name__] = typed_dict_shape(hint, special)
            return ref(hint.__name__)
        return None

    for root in (SynthesizeColorRulesRequest, SynthesizedColorRules):
        special(root)
    return {"types": types_}


def _indexed_color_rules(rules: Sequence[ColorRuleCondition]) -> dict:
    """The rules compiled as Generate compiles the color table, each colored by its index."""
    return build_specific_color_rules(DataFrame([
        dict(feature_type=r["feat"], qualifier_key=r["qual"], value=r["val"], color=str(i), caption="")
        for i, r in enumerate(rules)
    ]))


def _color_rule_ranks(item: RuleFeaturePayload, color_rules: dict) -> dict[int, int]:
    """The indexes of the rules a feature matches in precedence order (the first
    one colors it), each with its first rank."""
    feature = SimpleNamespace(type=item["type"], qualifiers=item["qualifiers"])
    ranks: dict[int, int] = {}
    for color, _, rank in iter_specific_color_rules(feature, color_rules, item.get("record", ""),
                                                    selector=item["selector"]):
        ranks.setdefault(int(color), rank)
    return ranks


def _winners(features: Sequence[RuleFeaturePayload], rules: Sequence[ColorRuleCondition]) -> list[int]:
    color_rules = _indexed_color_rules(rules)
    return [next(iter(_color_rule_ranks(item, color_rules)), -1) for item in features]


def evaluate_rules_json(features_json: str, rules_json: str, kind: str = "color") -> str:
    features = json.loads(features_json)
    rules = json.loads(rules_json)
    if kind == "color":
        color_rules = _indexed_color_rules(rules)
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
            ranks = _color_rule_ranks(item, color_rules)
            matched = list(ranks)
            priorities.append(list(ranks.values()))
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


_EXACT_SPECIAL = re.compile(r"[.*+?^${}()|[\]\\]")
# The qualifiers an exact condition reads besides a label text (Owner R-07).
_CONDITION_QUALIFIERS = ("product", "gene")
_REQUEST_DESCRIPTION = describe_rule_synthesis()

# An exact condition: (feature type, qualifier, pattern). A target of `hash`
# rules: (feature type, stable hash); duplicate record instances share one.
Condition = tuple[str, str, str]
HashTarget = tuple[str, str]


def exact_value_pattern(value: str) -> str:
    """The pattern that matches `value` exactly, escaped as the Web writes one
    (services/feature-selector.js `exactRegexValue`), so a rule either wrote
    reads as the same condition."""
    return "^" + _EXACT_SPECIAL.sub(lambda match: "\\" + match.group(0), value) + "$"


def _qualifier_rank(qualifier: str) -> int:
    order = [*_CONDITION_QUALIFIERS, *DEFAULT_LABEL_PRIORITY]
    return order.index(qualifier) if qualifier in order else len(order)


def _candidate_conditions(request: SynthesizeColorRulesRequest) -> list[Condition]:
    """The exact conditions of the selected features (Owner R-07): the feature
    type with its product, its gene, or a qualifier whose value is a label text
    the feature shows; product first, then gene, then the label qualifiers in
    Python's label priority, each in selection order."""
    found: dict[Condition, tuple[int, int]] = {}
    for position, (index, labels) in enumerate(zip(request["selected"], request["labels"])):
        item = request["features"][index]
        texts = {text.strip() for text in labels} - {""}
        for qualifier, values in item["qualifiers"].items():
            key = qualifier.lower()
            for value in values:
                if key != "hash" and value.strip() and (key in _CONDITION_QUALIFIERS or value.strip() in texts):
                    found.setdefault((item["type"], key, exact_value_pattern(value)), (_qualifier_rank(key), position))
    return sorted(found, key=found.__getitem__)


def _next_table(
    rules: Sequence[ColorRuleCondition],
    conditions: Sequence[ColorRuleRow],
    hash_rows: dict[HashTarget, ColorRuleRow],
    retired: dict[HashTarget, list[int]],
    anchors: dict[HashTarget, int],
) -> list[int | ColorRuleRow]:
    """The next table: the kept rows in order (by index); a condition before the
    first kept row of its type and qualifier, so it precedes them, its old rows
    of the same condition removed; a target's `hash` row in the place of its own
    old `hash` row (`retired`), else before the old `hash` row that matched it
    first (`anchors`), else last."""
    same_condition = {(row["feat"], row["qual"], row["val"]) for row in conditions}
    retired_at = {retired[key][0]: key for key in hash_rows if retired.get(key)}
    dropped = {position for rows in retired.values() for position in rows}
    placed_hashes = set(retired_at.values())
    table: list[int | ColorRuleRow] = []
    for position, rule in enumerate(rules):
        if (rule["feat"], rule["qual"].lower(), rule["val"]) in same_condition:
            continue
        if position in retired_at:
            table.append(hash_rows[retired_at[position]])
        elif position not in dropped:
            table.append(position)

    def insert_before(row: ColorRuleRow, kept: Callable[[int], bool]) -> None:
        at = next((at for at, entry in enumerate(table) if isinstance(entry, int) and kept(entry)), len(table))
        table.insert(at, row)

    for row in conditions:
        family = (row["feat"], row["qual"])
        insert_before(row, lambda index: (rules[index]["feat"], rules[index]["qual"].lower()) == family)
    for target, row in hash_rows.items():
        if target not in placed_hashes:
            anchor = anchors.get(target)
            insert_before(row, lambda index: index == anchor)
    return table


def synthesize_color_rules(request: SynthesizeColorRulesRequest) -> SynthesizedColorRules:
    """The fewest rules that give the selected features `color` and `caption`
    and leave every other feature's winning rule as it was (Owner R-07); with
    `scope` "features", one `hash` rule per feature.

    Python matches a single-feature rule by the stable hash, which duplicate
    record instances share: a feature that shares a selected feature's hash is a
    target too. A condition is taken only when every feature it matches is a
    target, largest first; the targets it does not cover take `hash` rules, and
    a target's own old `hash` rows give way. The table is then verified with the
    renderer's matcher: a target that a new row does not color (an older rule
    of an earlier qualifier wins) takes a `hash` row, and a table that would
    change the winner of any other feature is refused.
    """
    features, rules = request["features"], request["rules"]
    color, caption = request["color"], request["caption"]
    hash_of: list[HashTarget | None] = [
        (item["type"], item["selector"]["hash"]) if item["selector"]["hash"] else None for item in features
    ]
    chosen = {key for key in (hash_of[i] for i in request["selected"]) if key is not None}
    targets = set(request["selected"]) | {i for i, key in enumerate(hash_of) if key in chosen}

    old_rules = _indexed_color_rules(rules)
    old_matches = [list(_color_rule_ranks(item, old_rules)) for item in features]
    old_winners = [matched[0] if matched else -1 for matched in old_matches]

    candidates = _candidate_conditions(request) if request["scope"] == "group" else []
    candidate_rules = _indexed_color_rules([ColorRuleCondition(feat=t, qual=q, val=v) for t, q, v in candidates])
    hits: list[set[int]] = [set() for _ in candidates]
    for index, item in enumerate(features):
        for condition in _color_rule_ranks(item, candidate_rules):
            hits[condition].add(index)
    usable = sorted((c for c, found in enumerate(hits) if found and found <= targets), key=lambda c: -len(hits[c]))

    # The targets' hashes are the selected features' (`chosen`).
    spellings = {(kind, spelling): (kind, value) for kind, value in chosen
                 for spelling in (value, exact_value_pattern(value))}
    retired: dict[HashTarget, list[int]] = {}
    for position, rule in enumerate(rules):
        own = spellings.get((rule["feat"], rule["val"].strip())) if rule["qual"].lower() == "hash" else None
        if own:
            retired.setdefault(own, []).append(position)
    dropped = {position for rows in retired.values() for position in rows}

    forced: set[int] = set()
    for _attempt in range(3):
        covered: set[int] = set()
        conditions: list[ColorRuleRow] = []
        for c in usable:
            if hits[c] - covered - forced:
                feat, qual, val = candidates[c]
                conditions.append(ColorRuleRow(feat=feat, qual=qual, val=val, color=color, cap=caption))
                covered |= hits[c]
        needs_hash = {i for i in targets if i in forced or i not in covered}
        declined = {i for i in needs_hash if hash_of[i] is None}
        hash_rows: dict[HashTarget, ColorRuleRow] = {}
        anchors: dict[HashTarget, int] = {}
        for i in sorted(needs_hash - declined):
            key = hash_of[i]
            assert key is not None
            hash_rows.setdefault(key, ColorRuleRow(feat=key[0], qual="hash", val=key[1], color=color, cap=caption))
            anchor = next((r for r in old_matches[i] if r not in dropped and rules[r]["qual"].lower() == "hash"
                           and rules[r]["feat"] == key[0]), None)
            if anchor is not None:
                anchors.setdefault(key, anchor)
        table = _next_table(rules, conditions, hash_rows, retired, anchors)
        winners = _winners(features, [rules[row] if isinstance(row, int) else row for row in table])
        # The old row that colors each feature now; None for a new row, -1 for none.
        sources = [-1 if w < 0 else table[w] if isinstance(table[w], int) else None for w in winners]
        if any(sources[i] != old_winners[i] for i in range(len(features)) if i not in targets):
            raise ValueError("The rule synthesis would change the color of a feature that is not selected.")
        missed = {i for i in targets - declined if sources[i] is not None}
        if not missed:
            break
        if missed <= forced:
            raise ValueError("The rule synthesis could not color the selected features.")
        forced |= missed
    else:
        raise ValueError("The rule synthesis could not color the selected features.")
    return SynthesizedColorRules(
        table=table,
        exactRules=len([row for row in table if not isinstance(row, int) and row["qual"] != "hash"]),
        hashRules=len([row for row in table if not isinstance(row, int) and row["qual"] == "hash"]),
        declined=sorted(set(request["selected"]) & declined),
    )


def synthesize_color_rules_json(request_json: str) -> str:
    """:func:`synthesize_color_rules` for the Web's diagram helper `synthesizeColorRules`."""
    request = json.loads(request_json)
    # The feature payloads are read as `evaluate_rules_json` reads them; checking
    # each one's shape would cost as much as the synthesis on a large genome.
    features = request.get("features") if isinstance(request, dict) else None
    problems = violations({**request, "features": []} if isinstance(features, list) else request,
                          ref("SynthesizeColorRulesRequest"), _REQUEST_DESCRIPTION["types"], "request")
    if problems:
        raise ValueError("Unreadable rule synthesis request: " + "; ".join(problems[:3]))
    if len(request["labels"]) != len(request["selected"]) or not all(
            0 <= index < len(request["features"]) for index in request["selected"]):
        raise ValueError("Unreadable rule synthesis request: selected and labels do not match the features.")
    return json.dumps(synthesize_color_rules(request))
