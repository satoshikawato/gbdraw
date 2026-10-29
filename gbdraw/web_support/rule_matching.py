"""Evaluate browser rule drafts through the renderer's Python rule owners."""
from __future__ import annotations

import json
import re
from collections import defaultdict
from typing import Callable, Literal, Sequence, TypedDict, is_typeddict

from types import SimpleNamespace

from pandas import DataFrame

from gbdraw.contract_shapes import Shape, ref, typed_dict_shape, violations
from gbdraw.features.colors import build_specific_color_rules
from gbdraw.features.ids import feature_instance
from gbdraw.features.selector_values import ColorRuleList, iter_specific_color_rules
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
    """The drawn values a `hash`, `location`, or `record_location` rule reads,
    and the feature's catalog identity: `record_key` and `feature_id` (the
    catalog's ``recordKey`` and ``biologicalFeatureId``; None for a feature
    that is not a catalog feature). A `hash` rule reads the stable hash and
    the feature instance Python composes from that identity
    (``gbdraw.features.ids.feature_instance``)."""

    hash: str | None
    record_key: str | None
    feature_id: str | None
    location: str | None
    record_location: str | None


class RuleFeaturePayload(TypedDict):
    """A feature as Python's rule evaluation reads it (gbdraw/web/js/app/rule-matching.js `ruleFeaturePayload`)."""

    type: str
    record: str
    qualifiers: dict[str, list[str]]
    selector: RuleFeatureSelector
    label: str


class _MatchSelector(TypedDict):
    """The selector the renderer's matchers read (``iter_specific_color_rules``)."""

    hash: str | None
    instance: str | None
    location: str | None
    record_location: str | None


class _MatchFeature(TypedDict):
    """A payload feature with the selector the renderer's matchers read."""

    type: str
    record: str
    qualifiers: dict[str, list[str]]
    selector: _MatchSelector
    label: str


def _match_feature(item: RuleFeaturePayload) -> _MatchFeature:
    """The feature as Generate's matcher reads it: its instance composed from its
    catalog identity, None without a source hash or a catalog identity. Like the
    matchers, it reads a selector value the payload leaves out as None: the
    payloads are not checked one by one (`synthesize_color_rules_json`)."""
    selector = item["selector"]
    feature_hash, record_key, feature_id = selector.get("hash"), selector.get("record_key"), selector.get("feature_id")
    instance = feature_instance(feature_id, record_key) if feature_hash and record_key and feature_id else None
    return _MatchFeature(
        type=item["type"], record=item.get("record", ""), qualifiers=item["qualifiers"], label=item.get("label", ""),
        selector=_MatchSelector(hash=feature_hash, instance=instance, location=selector.get("location"),
                                record_location=selector.get("record_location")))


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
    Owner R-07); "features" colors each by itself (a `hash` rule naming the
    feature instance, as "This feature only" does); "reset" takes back the
    `hash` rows that name them and, when another rule still colors one and
    `color` is not empty, gives it `color` and `caption` by itself. `labels`
    holds, for each selected feature in the same order, the label texts it
    shows (its displayed and its annotation label). `rules` is the rule table
    at the write's turn, in order, with each row's color and caption.
    """

    scope: Literal["group", "features", "reset"]
    features: list[RuleFeaturePayload]
    selected: list[int]
    labels: list[list[str]]
    rules: list[ColorRuleRow]
    color: str
    caption: str


class RetargetedColorRule(TypedDict):
    """A copy of the request's row `row` with the pattern `val`."""

    row: int
    val: str


class SynthesizedColorRules(TypedDict):
    """The next rule table: a number keeps that row of the request's `rules`, a
    row with `row` copies that row under another pattern, any other object is a
    new row. A selected feature without a hash and without an exact condition
    takes no rule."""

    table: list[int | ColorRuleRow | RetargetedColorRule]


def describe_rule_synthesis() -> dict[str, dict[str, Shape]]:
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


IndexedColorRules = dict[str, dict[str, ColorRuleList]]


def _indexed_color_rules(rules: Sequence[ColorRuleCondition]) -> IndexedColorRules:
    """The rules compiled as Generate compiles the color table, each colored by its index."""
    return build_specific_color_rules(DataFrame([
        dict(feature_type=r["feat"], qualifier_key=r["qual"], value=r["val"], color=str(i), caption="")
        for i, r in enumerate(rules)
    ]))


def _color_rule_ranks(item: _MatchFeature, color_rules: IndexedColorRules) -> dict[int, int]:
    """The indexes of the rules a feature matches in precedence order (the first
    one colors it), each with its first rank."""
    feature = SimpleNamespace(type=item["type"], qualifiers=item["qualifiers"])
    ranks: dict[int, int] = {}
    selector = item["selector"]
    for color, _, rank in iter_specific_color_rules(feature, color_rules, item["record"],
                                                    selector=selector, instance=selector["instance"]):
        ranks.setdefault(int(color), rank)
    return ranks


def _winners(features: Sequence[_MatchFeature], rules: Sequence[ColorRuleCondition]) -> list[int]:
    color_rules = _indexed_color_rules(rules)
    return [next(iter(_color_rule_ranks(item, color_rules)), -1) for item in features]


def evaluate_rules_json(features_json: str, rules_json: str, kind: str = "color") -> str:
    features = [_match_feature(item) for item in json.loads(features_json)]
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
        selector = dict(item["selector"])
        record = item["record"]
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
                                             item["label"], label_rules, record, selector=selector)
            winners.append(int(winner) if winner is not None else -1)
    return json.dumps({"matches": matches, "winners": winners, "priorities": priorities})


_EXACT_SPECIAL = re.compile(r"[.*+?^${}()|[\]\\]")
# The qualifiers an exact condition reads besides a label text (Owner R-07).
_CONDITION_QUALIFIERS = ("product", "gene")
_REQUEST_DESCRIPTION = describe_rule_synthesis()

# An exact condition: (feature type, qualifier, pattern). A target of `hash`
# rules: (feature type, feature instance).
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
    in_place: dict[Condition, int],
) -> list[int | ColorRuleRow]:
    """The next table: the kept rows in order (by index); a condition in the
    place of the old row of the same condition named in `in_place`, else before
    the first kept row of its type and qualifier, so it precedes them, its old
    rows of the same condition removed; a target's `hash` row before the old `hash`
    row that matched it first (`anchors`) when that row comes before its own old
    `hash` row (a broader pattern would still color it), else in the place of
    its own old `hash` row (`retired`), else last."""
    same_condition = {(row["feat"], row["qual"], row["val"]) for row in conditions}
    placed_at = {in_place[key]: row for row in conditions
                 if (key := (row["feat"], row["qual"], row["val"])) in in_place}
    retired_at = {retired[key][0]: key for key in hash_rows
                  if retired.get(key) and anchors.get(key, len(rules)) > retired[key][0]}
    dropped = {position for rows in retired.values() for position in rows}
    placed_hashes = set(retired_at.values())
    table: list[int | ColorRuleRow] = []
    for position, rule in enumerate(rules):
        if (rule["feat"], rule["qual"].lower(), rule["val"]) in same_condition:
            if position in placed_at:
                table.append(placed_at[position])
            continue
        if position in retired_at:
            table.append(hash_rows[retired_at[position]])
        elif position not in dropped:
            table.append(position)

    def insert_before(row: ColorRuleRow, kept: Callable[[int], bool]) -> None:
        at = next((at for at, entry in enumerate(table) if isinstance(entry, int) and kept(entry)), len(table))
        table.insert(at, row)

    for row in conditions:
        if (row["feat"], row["qual"], row["val"]) in in_place:
            continue
        family = (row["feat"], row["qual"])
        insert_before(row, lambda index: (rules[index]["feat"], rules[index]["qual"].lower()) == family)
    for target, row in hash_rows.items():
        if target not in placed_hashes:
            anchor = anchors.get(target)
            insert_before(row, lambda index: index == anchor)
    return table


def _hash_targets(
    features: Sequence[_MatchFeature], selected: Sequence[int],
) -> tuple[list[HashTarget | None], set[int], dict[tuple[str, str], set[int]]]:
    """Each feature's `hash` target (its type and instance; None without an
    instance), the targets of an edit of `selected`, and the features of each
    type and stable hash. The list holds a feature once as drawn and once from
    its source, so the targets are every entry of a selected feature's target."""
    hash_of: list[HashTarget | None] = []
    holders: dict[tuple[str, str], set[int]] = defaultdict(set)
    for i, item in enumerate(features):
        instance = item["selector"]["instance"]
        hash_of.append((item["type"], instance) if instance else None)
        if item["selector"]["hash"]:
            holders[(item["type"], item["selector"]["hash"])].add(i)
    chosen = {key for key in (hash_of[i] for i in selected) if key is not None}
    return hash_of, set(selected) | {i for i, key in enumerate(hash_of) if key in chosen}, holders


def _spellings(value: str) -> tuple[str, str]:
    """A `hash` row value naming `value` exactly: as written, or as an exact pattern."""
    return value, exact_value_pattern(value)


def _reset_color_rules(request: SynthesizeColorRulesRequest, features: list[_MatchFeature]) -> SynthesizedColorRules:
    """The table without the `hash` rows that name a target: its instance rows,
    and its stable-hash rows. A stable-hash row that still colors other features
    (an older "This feature only" row on a feature drawn twice: the row wins for
    them) gives way to a copy per instance of each of those features, with its
    color and caption, so they and the Legend keep it; a feature that another
    row colors takes no copy. A target that another rule still colors
    then takes `color` and `caption` by itself, unless `color` is empty. A
    table that would change the winner of any other feature is refused."""
    rules = request["rules"]
    hash_of, targets, holders = _hash_targets(features, request["selected"])
    old_winners = _winners(features, rules)
    own: set[tuple[str, str]] = set()
    stable_of: dict[tuple[str, str], tuple[str, str]] = {}
    for i in targets:
        feat, stable = features[i]["type"], features[i]["selector"]["hash"]
        key = hash_of[i]
        own.update((feat, spelling) for spelling in (_spellings(key[1]) if key else ()))
        if stable:
            stable_of.update(((feat, spelling), (feat, stable)) for spelling in _spellings(stable))
    # Each entry of the next table, the row it reads as, and the old row it comes from.
    table: list[int | ColorRuleRow | RetargetedColorRule] = []
    rows: list[ColorRuleRow] = []
    origin: list[int] = []
    for position, rule in enumerate(rules):
        named = (rule["feat"], rule["val"].strip()) if rule["qual"].lower() == "hash" else None
        if named in own:
            continue
        if named in stable_of:
            others = {i for i in holders[stable_of[named]] - targets if old_winners[i] == position}
            if any(hash_of[i] is None for i in others):
                raise ValueError("The rule synthesis could not reset the color of the selected features.")
            for instance in sorted({key[1] for key in (hash_of[i] for i in others) if key}):
                table.append(RetargetedColorRule(row=position, val=exact_value_pattern(instance)))
                rows.append(ColorRuleRow(feat=rule["feat"], qual=rule["qual"], val=exact_value_pattern(instance),
                                         color=rule["color"], cap=rule["cap"]))
                origin.append(position)
            continue
        table.append(position)
        rows.append(rule)
        origin.append(position)
    winners = _winners(features, rows)
    if any((origin[winners[i]] if winners[i] >= 0 else -1) != old_winners[i]
           for i in range(len(features)) if i not in targets):
        raise ValueError("The rule synthesis would change the color of a feature that is not selected.")
    if not request["color"] or all(winners[i] < 0 for i in targets):
        return SynthesizedColorRules(table=table)
    colored = synthesize_color_rules(SynthesizeColorRulesRequest(
        scope="features", features=request["features"], selected=request["selected"],
        labels=[[] for _ in request["selected"]], rules=rows, color=request["color"], caption=request["caption"]))
    return SynthesizedColorRules(table=[table[entry] if isinstance(entry, int) else entry for entry in colored["table"]])


def synthesize_color_rules(request: SynthesizeColorRulesRequest) -> SynthesizedColorRules:
    """The fewest rules that give the selected features `color` and `caption`
    and leave every other feature's winning rule as it was (Owner R-07); with
    `scope` "features", one `hash` rule per feature; with `scope` "reset", the
    table without their own `hash` rows (`_reset_color_rules`).

    The targets are the selected features, with every other entry of the same
    feature instance (a feature is listed as drawn and from its source). A
    single-feature rule names the feature instance exactly, so a feature with
    the same stable hash (the same file drawn twice, identical features of one
    record) keeps its color. A condition is taken only when every feature it
    matches is a target, largest first; the targets it does not cover take
    `hash` rules, and a target's own old `hash` rows give way: its instance
    rows, and a stable-hash row that no other feature shares. The table is then
    verified with the renderer's matcher: a target that a new row does not
    color (an older rule of an earlier qualifier wins) takes a `hash` row, and
    a table that would change the winner of any other feature is refused.
    """
    features = [_match_feature(item) for item in request["features"]]
    if request["scope"] == "reset":
        return _reset_color_rules(request, features)
    rules = request["rules"]
    color, caption = request["color"], request["caption"]
    hash_of, targets, holders = _hash_targets(features, request["selected"])
    # The rows that name a target exactly: its instance, and its stable hash
    # when every feature with that hash is a target (a feature without an
    # instance gives way by its stable hash only).
    spellings: dict[tuple[str, str], HashTarget] = {}
    for i in sorted(targets):
        feat, key, stable = features[i]["type"], hash_of[i], features[i]["selector"]["hash"]
        names = [(key, key[1])] if key else []
        if stable and holders[(feat, stable)] <= targets:
            names.append((key or (feat, stable), stable))
        for target, value in names:
            for spelling in _spellings(value):
                spellings.setdefault((feat, spelling), target)

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
    # A condition that differs from the old row of the same condition at most
    # in its caption, where that row colored every feature the condition hits,
    # takes that row's place: the table order, and so the Legend order, stays
    # (OV-158; Owner 2026-10-29, OV-2189). The first such row is the one that
    # colors them.
    first_of: dict[Condition, int] = {}
    for position, rule in enumerate(rules):
        first_of.setdefault((rule["feat"], rule["qual"].lower(), rule["val"]), position)
    in_place = {candidates[c]: first_of[candidates[c]] for c in usable
                if candidates[c] in first_of and rules[first_of[candidates[c]]]["color"].lower() == color.lower()
                and all(old_winners[i] == first_of[candidates[c]] for i in hits[c])}
    # Such a condition comes first among those with as many hits, so another
    # condition that hits the same features (the product of a gene or label
    # row) does not take them and leave the old rule in place (OV-2191). The
    # sort is stable: without `in_place` the order is as before.
    usable.sort(key=lambda c: (-len(hits[c]), candidates[c] not in in_place))

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
            # An instance is named exactly: `f1@k` must not match `f1~1@k` or `f1@k:2`.
            hash_rows.setdefault(key, ColorRuleRow(feat=key[0], qual="hash", val=exact_value_pattern(key[1]),
                                                   color=color, cap=caption))
            anchor = next((r for r in old_matches[i] if r not in dropped and rules[r]["qual"].lower() == "hash"
                           and rules[r]["feat"] == key[0]), None)
            if anchor is not None:
                anchors.setdefault(key, anchor)
        table = _next_table(rules, conditions, hash_rows, retired, anchors, in_place)
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
    next_table: list[int | ColorRuleRow | RetargetedColorRule] = [*table]
    return SynthesizedColorRules(table=next_table)


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
