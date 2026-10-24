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

from gbdraw.features.selector_values import build_feature_selector_values
from gbdraw.web_support.rule_matching import evaluate_rules_json, exact_value_pattern, synthesize_color_rules_json

REPO = Path(__file__).resolve().parents[1]
COLOR = "#ff0000"


def _payload(feature: SeqFeature, record: str, label: str = "") -> dict:
    selector = build_feature_selector_values(feature, record)
    return dict(
        type=feature.type,
        record=record,
        qualifiers=selector["qualifiers"],
        selector=dict(hash=selector.get("hash"), location=selector.get("location"),
                      record_location=selector.get("record_location")),
        label=label,
    )


def _cds(index: int, **qualifiers: str) -> SeqFeature:
    return SeqFeature(SimpleLocation(index * 100, index * 100 + 60, strand=1), type="CDS",
                      qualifiers={key: [value] for key, value in qualifiers.items()})


def _synthesize(features: list[dict], selected: list[int], rules: list[dict] | None = None,
                labels: list[list[str]] | None = None, caption: str = "Group", scope: str = "group") -> dict:
    request = dict(scope=scope, features=features, selected=selected, labels=labels or [[] for _ in selected],
                   rules=[{key: rule[key] for key in ("feat", "qual", "val")} for rule in rules or []],
                   color=COLOR, caption=caption)
    return json.loads(synthesize_color_rules_json(json.dumps(request)))


def _next_table(answer: dict, current: list[dict]) -> list[dict]:
    return [current[row] if isinstance(row, int) else row for row in answer["table"]]


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
    assert answer["table"] == [0, {"feat": "CDS", "qual": "hash", "val": own, "color": COLOR, "cap": "Group"}]


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
    assert [row["val"] for row in table if row["qual"] == "hash"] == [features[4]["selector"]["hash"]]
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


def test_duplicate_record_instances_share_their_hash_and_are_one_target():
    feature = _cds(0, product="unique thing")
    other = _cds(1, product="unique thing")
    first, duplicate = _payload(feature, "rec"), _payload(feature, "rec")
    features = [first, duplicate, _payload(other, "rec")]
    answer = _synthesize(features, [0])
    table = _next_table(answer, [])
    assert [row["qual"] for row in table] == ["hash"]
    assert _colors(features, table) == [COLOR, COLOR, None]


def test_a_feature_without_hash_or_condition_is_declined():
    feature = _payload(_cds(0), "rec")
    feature["selector"]["hash"] = None
    answer = _synthesize([feature], [0])
    assert answer == {"table": []}


def test_a_malformed_request_is_refused():
    with pytest.raises(ValueError, match="selected"):
        synthesize_color_rules_json(json.dumps(dict(scope="group", features=[], labels=[], rules=[], color=COLOR, caption="x")))


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
    features[0]["selector"]["hash"], features[1]["selector"]["hash"] = "aaa1", "zzz9"
    current = [dict(feat="CDS", qual="hash", val="aaa1|zzz9", color="#0000ff", cap="both"),
               dict(feat="CDS", qual="hash", val="aaa1", color="#00ff00", cap="one")]
    answer = _synthesize(features, [0], rules=current, labels=[["kinase"]], scope=scope)
    table = _next_table(answer, current)
    assert _colors(features, table) == [COLOR, "#0000ff"]
    assert answer["table"] == [{"feat": "CDS", "qual": "hash", "val": "aaa1", "color": COLOR, "cap": "Group"}, 0]
