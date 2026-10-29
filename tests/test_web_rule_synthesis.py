"""OV-244: the Web's group color edits ask Python for the fewest rules that color
exactly the selected features (Owner R-07): exact rules on the feature type and
its product, gene, or label text that hit no unselected feature, and `hash`
rules for the rest."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation

from gbdraw.features.ids import feature_instance
from gbdraw.features.selector_values import build_feature_selector_values
from gbdraw.web_support.rule_matching import evaluate_rules_json, exact_value_pattern, synthesize_color_rules_json

REPO = Path(__file__).resolve().parents[1]
COLOR = "#ff0000"


def _payload(feature: SeqFeature, record: str, label: str = "", record_key: str | None = None,
             feature_id: str | None = None, catalog: bool = True) -> dict:
    """A feature as the Web sends it: `record_key` and `feature_id` are its
    catalog identity (the catalog's `recordKey` and `biologicalFeatureId`; by
    default `record` and the stable hash), which every catalog row carries;
    None without `catalog`, for a feature that is not a catalog feature."""
    selector = build_feature_selector_values(feature, record)
    feature_hash = selector.get("hash")
    return dict(
        type=feature.type,
        record=record,
        qualifiers=selector["qualifiers"],
        selector=dict(hash=feature_hash, record_key=(record_key or record) if catalog else None,
                      feature_id=(feature_id or feature_hash) if catalog else None,
                      location=selector.get("location"), record_location=selector.get("record_location")),
        label=label,
    )


def _instances(feature: SeqFeature, *record_keys: str, record: str = "rec") -> list[dict]:
    """The payloads of one feature drawn in each of `record_keys` (the same record twice)."""
    feature_hash = str(build_feature_selector_values(feature, record)["hash"])
    return [_payload(feature, record, record_key=key, feature_id=feature_hash) for key in record_keys]


def _instance(payload: dict) -> str:
    """The feature instance Python composes from a payload's catalog identity."""
    return feature_instance(payload["selector"]["feature_id"], payload["selector"]["record_key"])


def _cds(index: int, **qualifiers: str) -> SeqFeature:
    return SeqFeature(SimpleLocation(index * 100, index * 100 + 60, strand=1), type="CDS",
                      qualifiers={key: [value] for key, value in qualifiers.items()})


def _synthesize(features: list[dict], selected: list[int], rules: list[dict] | None = None,
                labels: list[list[str]] | None = None, caption: str = "Group", scope: str = "group",
                color: str = COLOR, joins: dict | None = None, drawn: list[int] | None = None) -> dict:
    request = dict(scope=scope, features=features, selected=selected, labels=labels or [[] for _ in selected],
                   rules=[{key: rule[key] for key in ("feat", "qual", "val", "color", "cap")} for rule in rules or []],
                   color=color, caption=caption, joins=joins, drawn=drawn)
    return json.loads(synthesize_color_rules_json(json.dumps(request)))


def _next_table(answer: dict, current: list[dict]) -> list[dict]:
    """The table as the Web writes it (rule-actions.js `groupColorRules`)."""
    return [current[row] if isinstance(row, int) else {**current[row["row"]], "val": row["val"]} if "row" in row
            else row for row in answer["table"]]


def _winners(features: list[dict], table: list[dict]) -> list[int]:
    return json.loads(evaluate_rules_json(json.dumps(features), json.dumps(table)))["winners"]


def _colors(features: list[dict], table: list[dict]) -> list[str | None]:
    return [table[w]["color"] if w >= 0 else None for w in _winners(features, table)]


def test_an_exact_pattern_is_written_as_the_web_writes_one():
    # The values of tests/web/feature-selector.test.mjs (`exactRegexValue`).
    assert exact_value_pattern("YP_009725295.1") == r"^YP_009725295\.1$"
    assert exact_value_pattern("  ATPase (A)+ [x]?  ") == r"^  ATPase \(A\)\+ \[x\]\?  $"


def test_a_whole_product_group_becomes_one_exact_rule():
    seqs = [_cds(0, product="hypothetical protein"), _cds(1, product="hypothetical protein"),
            _cds(2, product="DNA polymerase"), _cds(3, product="hypothetical protein")]
    features = [_payload(f, "rec") for f in seqs]
    answer = _synthesize(features, [0, 1, 3])
    assert answer["table"] == [{"feat": "CDS", "qual": "product", "val": "^hypothetical protein$",
                                "color": COLOR, "cap": "Group"}]


def test_a_feature_colored_by_itself_takes_a_hash_rule_in_place_of_its_own():
    seqs = [_cds(0, product="unique"), _cds(1, product="other")]
    features = [_payload(f, "rec") for f in seqs]
    own = features[0]["selector"]["hash"]
    current = [dict(feat="CDS", qual="product", val="^other$", color="#0000ff", cap="others"),
               dict(feat="CDS", qual="hash", val=own, color="#00ff00", cap="one")]
    answer = _synthesize(features, [0], rules=current, labels=[["unique"]], scope="features")
    assert answer["table"] == [0, {"feat": "CDS", "qual": "hash", "val": exact_value_pattern(_instance(features[0])),
                                   "color": COLOR, "cap": "Group"}]


def test_a_condition_that_hits_an_unselected_feature_is_refused():
    seqs = [_cds(0, product="hypothetical protein"), _cds(1, product="hypothetical protein"),
            _cds(2, product="Hypothetical protein")]  # rules match case-insensitively
    features = [_payload(f, "rec") for f in seqs]
    answer = _synthesize(features, [0, 1])
    table = _next_table(answer, [])
    assert [row["qual"] for row in table] == ["hash", "hash"]
    assert _colors(features, table) == [COLOR, COLOR, None]


def test_gene_and_label_text_conditions_cover_the_selection_and_the_rest_are_hash_rules():
    seqs = [_cds(0, product="kinase", gene="abcA"), _cds(1, product="kinase", gene="abcB"),
            _cds(2, product="kinase", gene="abcC", locus_tag="T_3"), _cds(3, product="other", locus_tag="T_4"),
            _cds(4, product="kinase"), _cds(5, product="other")]
    features = [_payload(f, "rec") for f in seqs]
    # Feature 3 shows its locus tag as its label; feature 4 shows a text no qualifier holds.
    answer = _synthesize(features, [0, 1, 3, 4], labels=[["kinase"], ["kinase"], ["T_4"], ["edited"]])
    table = _next_table(answer, [])
    conditions = sorted((row["qual"], row["val"]) for row in table if row["qual"] != "hash")
    assert conditions == [("gene", "^abcA$"), ("gene", "^abcB$"), ("locus_tag", "^T_4$")]
    assert [row["val"] for row in table if row["qual"] == "hash"] == [exact_value_pattern(_instance(features[4]))]
    assert _colors(features, table) == [COLOR, COLOR, None, COLOR, COLOR, None]


def test_no_unselected_feature_changes_its_winning_rule():
    seqs = [_cds(0, product="kinase", gene="k1"), _cds(1, product="kinase", gene="k2"),
            _cds(2, product="ligase", gene="l1")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="gene", val="^k", color="#00ff00", cap="k genes"),
               dict(feat="CDS", qual="product", val="ase$", color="#0000ff", cap="enzymes")]
    before = _winners(features, current)
    # The gene rule wins over a product rule for these features: the selected
    # feature takes a hash rule, and the others keep their rules.
    answer = _synthesize(features, [0], rules=current)
    table = _next_table(answer, current)
    winners = _winners(features, table)
    assert table[winners[0]]["color"] == COLOR
    assert [table[w] for w in winners[1:]] == [current[w] for w in before[1:]]
    assert answer["table"][:2] == [0, 1] and table[2]["qual"] == "hash"


def test_an_old_hash_rule_of_a_selected_feature_gives_way_to_its_new_rule():
    seqs = [_cds(0, product="kinase"), _cds(1, product="kinase"), _cds(2, product="ligase")]
    features = [_payload(f, "rec") for f in seqs]
    own = features[1]["selector"]["hash"]
    current = [dict(feat="CDS", qual="product", val="^ligase$", color="#0000ff", cap="ligases"),
               dict(feat="CDS", qual="hash", val=own, color="#00ff00", cap="one")]
    answer = _synthesize(features, [0, 1], rules=current)
    # A new condition goes before the first rule of its type and qualifier.
    assert answer["table"] == [{"feat": "CDS", "qual": "product", "val": "^kinase$", "color": COLOR, "cap": "Group"}, 0]


def test_a_rule_of_the_same_condition_gives_way_to_the_new_one_before_its_family():
    seqs = [_cds(0, product="kinase"), _cds(1, product="ligase")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="product", val="ase$", color="#0000ff", cap="enzymes"),
               dict(feat="CDS", qual="product", val="^kinase$", color="#00ff00", cap="kinases")]
    # In place, the broader rule before it would still color the feature.
    answer = _synthesize(features, [0], rules=current)
    assert answer["table"] == [{"feat": "CDS", "qual": "product", "val": "^kinase$", "color": COLOR, "cap": "Group"}, 0]


def test_a_rule_renamed_whole_keeps_its_place_in_the_table():
    # OV-2189 (Owner 2026-10-29): a new name for every feature of a rule's row
    # (the same condition and color, and the rule colored them) takes the
    # rule's place, so the table and Python's Legend keep their order (OV-158).
    seqs = [_cds(0, product="alpha protein"), _cds(1, product="dup alpha"), _cds(2, product="other")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha"),
               dict(feat="CDS", qual="product", val="^dup alpha$", color="#2a9d8f", cap="beta")]
    answer = _synthesize(features, [1], rules=current, color="#2A9D8F", caption="gamma")
    assert answer["table"] == [0, {"feat": "CDS", "qual": "product", "val": "^dup alpha$", "color": "#2A9D8F",
                                   "cap": "gamma"}]


@pytest.mark.parametrize(("name", "color", "extra"), [
    # A new color is a color edit, not a rename: the rule goes before its family, as before.
    ("a new color", "#000000", None),
    # The rule also colors a feature that is not selected: the condition is
    # refused and the selected feature takes a hash rule at the end; the other
    # feature keeps the old rule and its name.
    ("a rule that also colors an unselected feature", "#2a9d8f", "dup alpha"),
])
def test_only_a_rename_of_the_whole_rule_keeps_its_place(name, color, extra):
    seqs = [_cds(0, product="alpha protein"), _cds(1, product="dup alpha"), _cds(2, product=extra or "other")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha"),
               dict(feat="CDS", qual="product", val="^dup alpha$", color="#2a9d8f", cap="beta")]
    answer = _synthesize(features, [1], rules=current, color=color, caption="gamma")
    if extra is None:
        assert answer["table"] == [{"feat": "CDS", "qual": "product", "val": "^dup alpha$", "color": color,
                                     "cap": "gamma"}, 0], name
    else:
        table = _next_table(answer, current)
        assert answer["table"][:2] == [0, 1] and table[2]["qual"] == "hash", name
        assert [table[w]["cap"] for w in _winners(features, table)] == ["alpha", "gamma", "beta"], name


def test_a_rename_of_a_whole_locus_tag_row_appends_its_hash_rule_as_before():
    # A rule on a qualifier no condition reads (locus_tag, not a shown label),
    # whose only feature the edit renames: the feature takes a hash rule at the
    # end and the old rule stays in its place (Owner D-16 left this unchanged).
    seqs = [_cds(0, product="alpha protein"), _cds(1, product="x", locus_tag="T_1"), _cds(2, product="x", locus_tag="T_2")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="locus_tag", val="^T_1$", color="#2a9d8f", cap="beta"),
               dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha")]
    assert [current[w]["cap"] if w >= 0 else None for w in _winners(features, current)] == ["alpha", "beta", None]
    answer = _synthesize(features, [1], rules=current, color="#2a9d8f", caption="gamma")
    table = _next_table(answer, current)
    assert answer["table"][:2] == [0, 1] and len(table) == 3 and table[2]["qual"] == "hash"
    assert [table[w]["cap"] if w >= 0 else None for w in _winners(features, table)] == ["alpha", "gamma", None]


@pytest.mark.parametrize(("name", "beta", "beta_first", "products", "group", "row", "winners"), [
    # OV-2192 (Owner 2026-10-29, D-19 a): the Web says the edit takes the
    # features' Legend row whole (scope "row"); the new rows go at the place of
    # the rule that colored them, and the old rule stays for what it still colors.
    # A whole locus_tag row: the target's hash row.
    ("locus_tag row", ("locus_tag", "^T_1$"), True, ("x", "x"), [0, 1, "hash"], ["hash", 0, 1],
     ["alpha", "gamma", None]),
    # A whole locus_tag row whose feature has a product of its own: the new condition.
    ("locus_tag row, own product", ("locus_tag", "^T_1$"), False, ("unique", "x"), ["product", 0, 1],
     [0, "product", 1], ["alpha", "gamma", None]),
    # A product rule that also colors a feature the Web does not draw (T_2,
    # hidden): the Web takes the row whole; Python cannot tell it from a split
    # (scope "group"), where T_1 leaves the row.
    ("rule that also colors a hidden feature", ("product", "^dup alpha$"), False, ("dup alpha", "dup alpha"),
     [0, 1, "hash"], [0, "hash", 1], ["alpha", "gamma", "beta"]),
])
@pytest.mark.parametrize("scope", ["group", "row"])
def test_a_row_taken_whole_puts_its_new_rows_in_the_row_place(name, beta, beta_first, products, group, row, winners,
                                                              scope):
    alpha = dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha")
    beta_rule = dict(feat="CDS", qual=beta[0], val=beta[1], color="#2a9d8f", cap="beta")
    current = [beta_rule, alpha] if beta_first else [alpha, beta_rule]
    seqs = [_cds(0, product="alpha protein"), _cds(1, product=products[0], locus_tag="T_1"),
            _cds(2, product=products[1], locus_tag="T_2")]
    features = [_payload(f, "rec") for f in seqs]
    answer = _synthesize(features, [1], rules=current, color="#2a9d8f", caption="gamma", scope=scope)
    shape = [entry if isinstance(entry, int) else entry["qual"] for entry in answer["table"]]
    assert shape == (group if scope == "group" else row), name
    # The old rule stays in the table for every feature it still colors.
    table = _next_table(answer, current)
    assert beta_rule in table, name
    assert [table[w]["cap"] if w >= 0 else None for w in _winners(features, table)] == winners, name


@pytest.mark.parametrize("scope", ["group", "row"])
def test_a_same_pattern_rename_is_the_same_for_a_row_taken_whole(scope):
    # The same-pattern re-caption (OV-2189) does not depend on the scope.
    seqs = [_cds(0, product="alpha protein"), _cds(1, gene="dupA", product="dup alpha"), _cds(2, product="other")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha"),
               dict(feat="CDS", qual="gene", val="^dupA$", color="#2a9d8f", cap="beta")]
    answer = _synthesize(features, [1], rules=current, color="#2a9d8f", caption="gamma", scope=scope)
    assert answer["table"] == [0, {**current[1], "cap": "gamma"}]


@pytest.mark.parametrize("beta_first", [True, False], ids=["beta-row-first", "alpha-row-first"])
@pytest.mark.parametrize(("qualifier", "value", "labels", "gene_first"), [
    # OV-2191: a gene rule row; the feature also has a product, which hits the
    # same features and used to be taken first, so the old rule stayed.
    ("gene", "dupA", [], True),
    ("gene", "dupA", [], False),
    # A label qualifier the feature shows as its label.
    ("note", "beta note", ["beta note"], True),
], ids=["gene-before-product", "product-before-gene", "note-label"])
def test_a_gene_or_label_rule_renamed_whole_keeps_its_place_in_the_table(qualifier, value, labels, gene_first,
                                                                          beta_first):
    # The same-pattern rule is re-captioned in place even when the feature's
    # product would hit the same features (Owner D-13 / D-16, OV-2189).
    named = {qualifier: value, "product": "dup alpha"} if gene_first else {"product": "dup alpha", qualifier: value}
    seqs = [_cds(0, product="alpha protein"), _cds(1, **named), _cds(2, product="other")]
    features = [_payload(f, "rec") for f in seqs]
    beta = dict(feat="CDS", qual=qualifier, val=exact_value_pattern(value), color="#2a9d8f", cap="beta")
    alpha = dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha")
    current = [beta, alpha] if beta_first else [alpha, beta]
    answer = _synthesize(features, [1], rules=current, labels=[labels], color="#2a9d8f", caption="gamma")
    renamed = {**beta, "cap": "gamma"}
    assert answer["table"] == ([renamed, 1] if beta_first else [0, renamed])


@pytest.mark.parametrize("scope", ["features", "group"])
def test_one_of_two_record_instances_with_the_same_hash_is_the_only_target(scope):
    # PD-OI-069 revision 2: the same file drawn twice; its twin keeps its color.
    first, twin = _instances(_cds(0, product="unique thing"), "linear-seq-a", "linear-seq-b")
    features = [first, twin, _payload(_cds(1, product="other"), "rec")]
    answer = _synthesize(features, [0], labels=[["unique thing"]], scope=scope)
    table = _next_table(answer, [])
    assert table == [{"feat": "CDS", "qual": "hash", "val": exact_value_pattern(_instance(first)),
                      "color": COLOR, "cap": "Group"}]
    assert _colors(features, table) == [COLOR, None, None]


def test_one_of_identical_features_of_a_record_is_the_only_target():
    feature = _cds(0)
    feature_hash = build_feature_selector_values(feature, "rec")["hash"]
    features = [_payload(feature, "rec", record_key="k", feature_id=f"{feature_hash}~{n}") for n in (1, 2)]
    table = _next_table(_synthesize(features, [1], scope="features"), [])
    assert _colors(features, table) == [None, COLOR]


def test_a_hash_rule_a_twin_shares_stays_and_the_new_instance_row_wins_before_it():
    first, twin = _instances(_cds(0), "k1", "k2")
    features = [first, twin]
    current = [dict(feat="CDS", qual="hash", val=first["selector"]["hash"], color="#0000ff", cap="both")]
    answer = _synthesize(features, [0], rules=current, scope="features")
    assert answer["table"] == [{"feat": "CDS", "qual": "hash", "val": exact_value_pattern(_instance(first)),
                                "color": COLOR, "cap": "Group"}, 0]
    assert _colors(features, _next_table(answer, current)) == [COLOR, "#0000ff"]


def test_a_hash_rule_of_a_feature_without_a_twin_gives_way_to_its_instance_row():
    (only,) = _instances(_cds(0), "k1")
    current = [dict(feat="CDS", qual="hash", val=only["selector"]["hash"], color="#0000ff", cap="old")]
    answer = _synthesize([only], [0], rules=current, scope="features")
    assert answer["table"] == [{"feat": "CDS", "qual": "hash", "val": exact_value_pattern(_instance(only)),
                                "color": COLOR, "cap": "Group"}]


def test_a_feature_listed_as_drawn_and_from_its_source_is_one_target():
    # The Web lists each feature once as drawn and once from its source; both
    # entries are the same instance, and the twin of the other record keeps its color.
    first, twin = _instances(_cds(0, product="unique thing"), "k1", "k2")
    features = [first, twin, dict(first), dict(twin)]
    table = _next_table(_synthesize(features, [0], scope="features"), [])
    assert _colors(features, table) == [COLOR, None, COLOR, None]


def test_a_feature_without_a_catalog_identity_takes_no_hash_rule():
    features = [_payload(_cds(0), "rec", catalog=False)]
    assert _synthesize(features, [0], scope="features")["table"] == []


def test_a_reset_of_one_copy_retargets_the_hash_rule_its_twin_shares():
    # An older "This feature only" row names the stable hash both copies share;
    # resetting one copy keeps the twin's color and Legend row.
    first, twin = _instances(_cds(0), "k1", "k2")
    features = [first, twin, dict(first), dict(twin)]
    shared = dict(feat="CDS", qual="hash", val=first["selector"]["hash"], color="#0000ff", cap="both", fromFile=True)
    current = [shared, dict(feat="CDS", qual="product", val="^x$", color="#00ff00", cap="x")]
    answer = _synthesize(features, [0], rules=current, scope="reset", color="")
    assert answer["table"] == [{"row": 0, "val": exact_value_pattern(_instance(twin))}, 1]
    table = _next_table(answer, current)
    assert table[0] == {**shared, "val": exact_value_pattern(_instance(twin))}
    assert _colors(features, table) == [None, "#0000ff", None, "#0000ff"]


def test_a_reset_of_one_copy_narrows_the_shared_hash_rule_to_the_copies_it_still_colors():
    # The twin wins by its own instance row, so the shared row colors only the
    # reset copy; the reset leaves no row behind that colors nothing.
    first, twin = _instances(_cds(0), "k1", "k2")
    current = [dict(feat="CDS", qual="hash", val=exact_value_pattern(_instance(first)), color="#ff0000", cap="red"),
               dict(feat="CDS", qual="hash", val=first["selector"]["hash"], color="#999999", cap="grey")]
    answer = _synthesize([first, twin], [1], rules=current, scope="reset", color="")
    assert answer["table"] == [0]
    assert _colors([first, twin], _next_table(answer, current)) == ["#ff0000", None]


def test_a_reset_takes_back_the_instance_row_and_a_hash_row_no_other_feature_holds():
    first, twin = _instances(_cds(0), "k1", "k2")
    (other,) = _instances(_cds(1), "k1")
    current = [dict(feat="CDS", qual="hash", val=exact_value_pattern(_instance(first)), color="#0000ff", cap="mine"),
               dict(feat="CDS", qual="hash", val=other["selector"]["hash"], color="#00ff00", cap="other"),
               dict(feat="CDS", qual="hash", val=exact_value_pattern(_instance(twin)), color="#0000ff", cap="mine")]
    assert _synthesize([first, twin, other], [0], rules=current, scope="reset", color="")["table"] == [1, 2]
    assert _synthesize([first, twin, other], [2], rules=current, scope="reset", color="")["table"] == [0, 2]


def test_a_reset_feature_another_rule_still_colors_takes_the_reset_color():
    first, twin = _instances(_cds(0, product="x"), "k1", "k2")
    current = [dict(feat="CDS", qual="hash", val=first["selector"]["hash"], color="#0000ff", cap="both"),
               dict(feat="CDS", qual="product", val="^x$", color="#00ff00", cap="x")]
    answer = _synthesize([first, twin], [0], rules=current, scope="reset", color="#999999", caption="other proteins")
    table = _next_table(answer, current)
    assert _colors([first, twin], table) == ["#999999", "#0000ff"]
    assert {"feat": "CDS", "qual": "hash", "val": exact_value_pattern(_instance(first)), "color": "#999999",
            "cap": "other proteins"} in table
    # Without another rule, the feature takes no row: it shows its palette color.
    answer = _synthesize([first, twin], [0], rules=current[:1], scope="reset", color="#999999", caption="other proteins")
    assert _colors([first, twin], _next_table(answer, current[:1])) == [None, "#0000ff"]


def test_a_feature_without_hash_or_condition_is_declined():
    feature = _payload(_cds(0), "rec")
    feature["selector"]["hash"] = None
    answer = _synthesize([feature], [0])
    assert answer == {"table": [], "legend": None}


def test_a_malformed_request_is_refused():
    with pytest.raises(ValueError, match="selected"):
        synthesize_color_rules_json(json.dumps(dict(scope="group", features=[], labels=[], rules=[], color=COLOR, caption="x",
                                                    joins=None, drawn=None)))
    with pytest.raises(ValueError, match="joined row"):
        synthesize_color_rules_json(json.dumps(dict(scope="group", features=[], selected=[], labels=[], rules=[],
                                                    color=COLOR, caption="x", joins=dict(at=0, rows=[]), drawn=None)))
    with pytest.raises(ValueError, match="drawn feature"):
        synthesize_color_rules_json(json.dumps(dict(scope="group", features=[], selected=[], labels=[], rules=[],
                                                    color=COLOR, caption="x", joins=None, drawn=[0])))


def _genbank_payloads(path: Path) -> list[dict]:
    return [_payload(f, record.id) for record in SeqIO.parse(path, "genbank")
            for f in record.features if f.type == "CDS"]


@pytest.mark.parametrize(("path", "count"), [
    (REPO / "tests" / "test_inputs" / "AP027078.gb", 236),
    (REPO / "examples" / "NC_010162.gb", 2761),
])
def test_every_hypothetical_protein_of_a_genome_becomes_one_rule(path, count):
    features = _genbank_payloads(path)
    products = [f["qualifiers"].get("product", [""])[0] for f in features]
    selected = [i for i, product in enumerate(products) if product == "hypothetical protein"]
    assert len(selected) == count
    assert Counter(p.lower() for p in products)["hypothetical protein"] == count
    answer = _synthesize(features, selected, labels=[["hypothetical protein"] for _ in selected],
                         caption="hypothetical protein")
    assert answer["table"] == [{"feat": "CDS", "qual": "product", "val": "^hypothetical protein$",
                                "color": COLOR, "cap": "hypothetical protein"}]


@pytest.mark.parametrize("scope", ["features", "group"])
def test_a_new_hash_rule_goes_before_a_broader_old_hash_rule_that_colors_the_feature(scope):
    # Review case: a broader `hash` row before the feature's own `hash` row
    # colors it, so the new row cannot take the own row's place.
    seqs = [_cds(0, product="kinase"), _cds(1, product="kinase")]
    features = [_payload(f, "rec") for f in seqs]
    for feature, value in zip(features, ("aaa1", "zzz9")):
        feature["selector"]["hash"] = feature["selector"]["feature_id"] = value
    current = [dict(feat="CDS", qual="hash", val="aaa1|zzz9", color="#0000ff", cap="both"),
               dict(feat="CDS", qual="hash", val="aaa1", color="#00ff00", cap="one")]
    answer = _synthesize(features, [0], rules=current, labels=[["kinase"]], scope=scope)
    table = _next_table(answer, current)
    assert _colors(features, table) == [COLOR, "#0000ff"]
    assert answer["table"] == [{"feat": "CDS", "qual": "hash", "val": r"^aaa1@rec$", "color": COLOR, "cap": "Group"}, 0]


TEAL = "#2a9d8f"


def _legend_order(features: list[dict], table: list[dict]) -> list[str]:
    """The rule rows Python's Legend draws, in its order: each caption and color
    a feature takes, at its first rule in the table (prepare_legend_table)."""
    used = {(table[w]["cap"], table[w]["color"]) for w in _winners(features, table) if w >= 0}
    return list(dict.fromkeys(rule["cap"] for rule in table if (rule["cap"], rule["color"]) in used))


def _join_seqs() -> list[SeqFeature]:
    return [_cds(0, product="alpha protein"), _cds(1, product="dup alpha"), _cds(2, product="beta protein"),
            _cds(3, product="dup beta")]


def _join_case(beta: str) -> tuple[list[dict], list[dict]]:
    """A beta row of pattern `beta`, a delta row, and an alpha row of beta's color."""
    seqs = _join_seqs()
    current = [dict(feat="CDS", qual="product", val=beta, color=TEAL, cap="beta"),
               dict(feat="CDS", qual="product", val="^beta protein$", color="#e9c46a", cap="delta"),
               dict(feat="CDS", qual="product", val="^alpha protein$", color=TEAL, cap="alpha")]
    return [_payload(f, "rec") for f in seqs], current


def test_a_rename_that_joins_a_row_through_a_same_pattern_rule_keeps_the_joined_row_in_place():
    # OV-2193 (Owner 2026-10-29): the re-captioned rule goes just before the
    # joined row's first rule, so Python's Legend keeps the joined row where the
    # live Legend keeps it; without `joins` it stays in its own place (OV-2189).
    features, current = _join_case("^dup alpha$")
    joined = dict(feat="CDS", qual="product", val="^dup alpha$", color=TEAL, cap="alpha")
    alone = _synthesize(features, [1], rules=current, color=TEAL, caption="alpha")
    assert alone["table"] == [joined, 1, 2]
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="alpha", joins=dict(at=2, rows=[]))
    assert answer["table"] == [1, joined, 2]
    assert _legend_order(features, _next_table(answer, current)) == ["delta", "alpha"]


def test_a_join_of_part_of_a_row_takes_a_hash_rule_before_the_joined_row():
    # The feature's own condition would lose to the broader beta rule before the
    # joined row's first rule, so it takes a `hash` rule there instead.
    features, current = _join_case("^dup")
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="alpha", joins=dict(at=2, rows=[]))
    assert [entry if isinstance(entry, int) else entry["qual"] for entry in answer["table"]] == [0, 1, "hash", 2]
    table = _next_table(answer, current)
    assert [table[w]["cap"] for w in _winners(features, table)] == ["alpha", "alpha", "delta", "beta"]
    assert _legend_order(features, table) == ["beta", "delta", "alpha"]


def test_a_row_of_own_hash_rules_that_joins_a_row_goes_before_the_joined_row():
    # A "This feature only" row: its own `hash` row gives way, and the new rule
    # goes before the joined row's first rule, not in the old row's place.
    features, current = _join_case("^dup alpha$")
    current[0] = dict(feat="CDS", qual="hash", val=exact_value_pattern(_instance(features[1])), color=TEAL, cap="beta")
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="alpha", joins=dict(at=2, rows=[]))
    assert answer["table"] == [1, dict(feat="CDS", qual="product", val="^dup alpha$", color=TEAL, cap="alpha"), 2]


def test_rules_renamed_onto_a_listed_row_move_just_before_its_first_rule():
    # The Legend editor's rename: the Web re-captions the row's rules, and
    # Python moves the ones before the joined row's first rule to just before it.
    features, current = _join_case("^dup")
    renamed = [{**current[0], "cap": "alpha"}, *current[1:]]
    answer = _synthesize(features, [], rules=renamed, color=TEAL, caption="alpha", joins=dict(at=2, rows=[0]))
    assert answer["table"] == [1, 0, 2]
    assert _legend_order(features, _next_table(answer, renamed)) == ["delta", "alpha"]


def test_rules_renamed_onto_a_listed_row_stay_in_place_when_the_move_would_recolor_a_feature():
    # Moved past the delta rule, the renamed rule would lose `dup beta` to it:
    # the rows keep their places, as without a join (the joined row moves).
    features, current = _join_case("^dup")
    renamed = [{**current[0], "cap": "alpha"}, {**current[1], "val": "beta"}, current[2]]
    answer = _synthesize(features, [], rules=renamed, color=TEAL, caption="alpha", joins=dict(at=2, rows=[0]))
    assert answer["table"] == [0, 1, 2]


def test_a_same_pattern_rule_after_the_joined_row_stays_in_place():
    # Review M1 of #1467: re-captioned in its own place, the rule already keeps
    # the joined row's place (it stands after the joined row's first rule), so
    # the join does not move it (OV-2189).
    seqs = [_cds(0, product="alpha protein"), _cds(1, product="dup alpha"), _cds(2, product="beta protein"),
            _cds(3, product="dup beta")]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="product", val="^alpha protein$", color=TEAL, cap="alpha"),
               dict(feat="CDS", qual="product", val="^beta protein$", color="#e9c46a", cap="delta"),
               dict(feat="CDS", qual="product", val="^dup alpha$", color=TEAL, cap="beta"),
               dict(feat="CDS", qual="product", val="^dup beta$", color=TEAL, cap="beta")]
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="alpha", joins=dict(at=0, rows=[]))
    assert answer["table"] == [0, 1, {**current[2], "cap": "alpha"}, 3]


def test_an_own_hash_row_after_the_joined_row_keeps_its_place():
    # The same for a "This feature only" row after the joined row: the new
    # `hash` row takes the place of the old one.
    seqs = [_cds(0, product="alpha protein"), _cds(1, product="dup"), _cds(2, product="dup")]
    features = [_payload(f, "rec") for f in seqs]
    own = dict(feat="CDS", qual="hash", val=exact_value_pattern(_instance(features[1])), color=TEAL, cap="beta")
    current = [dict(feat="CDS", qual="product", val="^alpha protein$", color=TEAL, cap="alpha"),
               dict(feat="CDS", qual="product", val="^beta$", color="#e9c46a", cap="delta"), own]
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="alpha", joins=dict(at=0, rows=[]))
    assert answer["table"] == [0, 1, {**own, "cap": "alpha"}]


def _generate_legend(seqs: list[SeqFeature], table: list[dict]) -> list[str]:
    """The feature rows Generate draws from `table` for a record of `seqs`, by
    name: Generate's presence check (`feature_legend_presence`) on the record,
    then its Legend owner (`prepare_legend_table`)."""
    from types import SimpleNamespace

    import pandas as pd
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.features.colors import build_specific_color_rules
    from gbdraw.legend.sources import FeatureLegendSource, feature_legend_presence
    from gbdraw.legend.table import prepare_legend_table

    record = SeqRecord(Seq("A" * (max(int(f.location.end) for f in seqs) + 100)), id="rec", features=list(seqs))
    color_table = pd.DataFrame([[r["feat"], r["qual"], r["val"], r["color"], r["cap"]] for r in table],
                               columns=["feature_type", "qualifier_key", "value", "color", "caption"]) if table else None
    source = FeatureLegendSource(
        color_table=color_table,
        default_colors=pd.DataFrame([["default", "#d3d3d3"]], columns=["feature_type", "color"]),
        block_stroke_color="none", block_stroke_width=0,
        selected_features=tuple(dict.fromkeys(f.type for f in seqs)), feature_visibility_rules=None,
        specific_color_rules=build_specific_color_rules(color_table), default_color_map={"default": "#d3d3d3"},
        record_features=())
    present, used, default_used = feature_legend_presence([record], source)
    numeric = SimpleNamespace(dinucleotide="GC", high_fill_color="#000000", low_fill_color="#000000",
                              stroke_color="none", stroke_width=0)
    legend = prepare_legend_table(numeric, numeric, source, present, used_color_rules=used,
                                  default_used_features=default_used, show_gc=False, show_skew=False, show_depth=False)
    return [str(names.name) for names in legend.names().values()]


def _answer_legend(answer: dict, table: list[dict]) -> list[str]:
    """The answer's Legend order by name: a rule row by its rule's caption."""
    return [table[row]["cap"] if isinstance(row, int) else row for row in answer["legend"]]


def test_the_answer_gives_the_legend_order_of_the_next_table_where_the_old_row_moves():
    # OV-2247: a split whose old row also moves in Python's Legend. The feature's
    # rule is re-captioned in place, first in the table, so the old row now
    # starts at its later rule, after alpha. The answer gives that order, so
    # the live Legend shows it at once.
    seqs = _join_seqs()
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="CDS", qual="product", val="^dup alpha$", color=TEAL, cap="beta"),
               dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha"),
               dict(feat="CDS", qual="product", val="^dup beta$", color=TEAL, cap="beta")]
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="gamma", drawn=[0, 1, 2, 3])
    table = _next_table(answer, current)
    assert [rule["cap"] for rule in table] == ["gamma", "alpha", "beta"]
    assert answer["legend"] == [0, 1, 2, "other proteins"]
    assert _answer_legend(answer, table) == _generate_legend(seqs, table)


def test_the_answer_gives_the_old_order_when_the_join_is_refused():
    # OV-2248 (Owner 2026-10-29: Generate keeps the old order): moving the
    # renamed rule past the delta rule would recolor `dup beta`, so the table
    # stays and the renamed row keeps its first rule's place, before delta.
    features, current = _join_case("^dup")
    renamed = [{**current[0], "cap": "alpha"}, {**current[1], "val": "beta"}, current[2]]
    answer = _synthesize(features, [], rules=renamed, color=TEAL, caption="alpha", joins=dict(at=2, rows=[0]),
                         drawn=[0, 1, 2, 3])
    assert answer["table"] == [0, 1, 2]
    assert answer["legend"] == [0, 1]
    assert _answer_legend(answer, renamed) == _generate_legend(_join_seqs(), renamed) == ["alpha", "delta"]


def test_the_answer_groups_the_legend_rows_by_feature_type_in_drawing_order():
    # Review M2 of #1467: Python draws the rule rows of the first drawn type
    # first, whatever their table order: a tRNA split off its row goes after
    # the CDS rows when a CDS is drawn first.
    trna = SeqFeature(SimpleLocation(900, 970, strand=1), type="tRNA", qualifiers={"product": ["tRNA-Gly"]})
    seqs = [_cds(0, product="alpha protein"), _cds(1, product="dup beta"), trna]
    features = [_payload(f, "rec") for f in seqs]
    current = [dict(feat="tRNA", qual="product", val="^tRNA-Gly$", color=TEAL, cap="beta"),
               dict(feat="CDS", qual="product", val="^dup beta$", color=TEAL, cap="beta"),
               dict(feat="CDS", qual="product", val="^alpha protein$", color="#e63946", cap="alpha")]
    answer = _synthesize(features, [2], rules=current, color=TEAL, caption="gamma", drawn=[0, 1, 2])
    table = _next_table(answer, current)
    assert answer["legend"] == [1, 2, 0]
    assert _answer_legend(answer, table) == _generate_legend(seqs, table) == ["beta", "alpha", "gamma"]


def test_the_answer_reads_only_the_drawn_features():
    # A feature the Result does not draw (another Result's, or a source copy)
    # takes no row; a drawn type no rule colors keeps its own row.
    features, current = _join_case("^dup alpha$")
    answer = _synthesize(features, [1], rules=current, color=TEAL, caption="gamma", drawn=[1, 3])
    table = _next_table(answer, current)
    assert _answer_legend(answer, table) == _generate_legend([_join_seqs()[1], _join_seqs()[3]], table) == [
        "gamma", "other proteins"]
    assert _synthesize(features, [1], rules=current, color=TEAL, caption="gamma")["legend"] is None
