"""OV-193: a long rule list is looked up through its literal index, with the
results of testing every rule."""
import random
import re
import string
from types import SimpleNamespace

from pandas import DataFrame

from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.selector_values import ColorRuleList, iter_specific_color_rules

_ALPHABET = "abcAB.-$^*_ "
_VALUES = ["", "\n", "abc\n", "İb", "ſa", "Ka", "a.b", "A-B", "x$", "β"]


def _js_escape(text):
    return re.sub(r"[.*+?^${}()|[\]\\]", lambda match: "\\" + match.group(0), text)


def _pattern(rng):
    literal = "".join(rng.choice(_ALPHABET) for _ in range(rng.randint(0, 3)))
    escape = rng.choice([re.escape, _js_escape])
    return rng.choice([
        escape(literal), f"^{escape(literal)}", f"{escape(literal)}$", f"^{escape(literal)}$",
        literal.swapcase().replace("$", "\\$"), "a.c", "[ab]+", "(?i)ab", "a|b", "", "^", "$", "^$",
        "\\d", "ab\\$", "a\\\\$", "a{1}", "ſ", "K",
    ])


def _value(rng):
    if rng.random() < 0.15:
        return rng.choice(_VALUES)
    return "".join(rng.choice(_ALPHABET) for _ in range(rng.randint(0, 6)))


def _plain(color_map):
    return {feature_type: {key: list(rules) for key, rules in by_key.items()}
            for feature_type, by_key in color_map.items()}


def test_indexed_lookup_yields_what_testing_every_rule_yields():
    rng = random.Random(193)
    rows = []
    for _ in range(400):
        pattern = _pattern(rng)
        try:
            re.compile(pattern)
        except re.error:
            continue
        rows.append(dict(feature_type=rng.choice(["CDS", "CDS", "*", "tRNA"]),
                         qualifier_key=rng.choice(["hash", "product", "location", "record_location"]),
                         value=pattern, color=str(len(rows)), caption=""))
    rows += [dict(rows[i], color=str(len(rows) + i)) for i in range(20)]  # duplicates
    color_map, _ = preprocess_color_tables(DataFrame(rows), DataFrame(columns=["feature_type", "color"]))
    assert any(rules.literal_index is not None for by_key in color_map.values() for rules in by_key.values())
    plain = _plain(color_map)
    for _ in range(600):
        feature = SimpleNamespace(type=rng.choice(["CDS", "tRNA", "rRNA"]),
                                  qualifiers={"product": [_value(rng) for _ in range(rng.randint(0, 3))]})
        selector = {key: _value(rng) for key in ["hash", "location", "record_location"]}
        assert list(iter_specific_color_rules(feature, color_map, "rec", selector=selector)) == list(
            iter_specific_color_rules(feature, plain, "rec", selector=selector))


class _CountingPattern:
    def __init__(self, source, counter):
        self._compiled = re.compile(source, re.IGNORECASE)
        self.pattern, self.flags, self._counter = source, self._compiled.flags, counter

    def search(self, value):
        self._counter[0] += 1
        return self._compiled.search(value)


def test_hash_rules_test_only_the_rules_a_feature_can_match():
    for count in (50, 200):
        searches = [0]
        hashes = [f"{index:016x}" for index in range(count)]
        color_map = {"CDS": {"hash": ColorRuleList(
            [(_CountingPattern(value, searches), str(index), "") for index, value in enumerate(hashes)]
            + [(_CountingPattern(f"^{value}$", searches), str(count + index), "")
               for index, value in enumerate(hashes)]
        )}}
        for value in hashes:
            feature = SimpleNamespace(type="CDS", qualifiers={})
            colors = [color for color, _, _ in iter_specific_color_rules(
                feature, color_map, selector={"hash": value.upper()})]
            assert colors == [str(hashes.index(value)), str(count + hashes.index(value))]
        assert searches[0] == 2 * count


def test_a_short_list_and_a_non_ascii_value_test_every_rule():
    assert ColorRuleList([(re.compile("a"), "1", "")] * 31).literal_index is None
    searches = [0]
    rules = ColorRuleList([(_CountingPattern(f"^{letter}$", searches), letter, "")
                           for letter in string.ascii_letters[:40]])
    feature = SimpleNamespace(type="CDS", qualifiers={})
    assert [color for color, _, _ in iter_specific_color_rules(
        feature, {"CDS": {"hash": rules}}, selector={"hash": "K"})] == ["k", "K"]
    assert searches[0] == 40
