"""Per-feature edits addressed by original-source identity (design Q4, PR-Q4-2)."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw import (
    FeatureOptions,
    FeaturePlacementOverride,
    FeaturePlacementTarget,
    LinearOptions,
    draw_linear,
    read_genbank,
)
from gbdraw.annotations import AnnotationOptions, AnnotationSet, FeatureIdentitySpan, RegionAnnotation
from gbdraw.api.options import CircularDiagramOptions, CircularMultiRecordOptions, LinearDiagramOptions
from gbdraw.api.request_render import (
    _extract_linear_request_proteins,
    build_request_diagram,
    plan_request,
)
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
)
from gbdraw.features.overrides import FeatureOverride
from gbdraw.features.selector_values import get_feature_hash
from gbdraw.features.source import FeatureIdentity, resolve_feature_identities
from gbdraw.io.regions import parse_region_spec

FIXTURES = Path(__file__).parent / "fixtures"
COLLIDE = FIXTURES / "b_collide.gb"
TYPES = ("CDS", "tRNA", "misc_feature")
# b_collide.gb draws feature Y of the 201-4000 crop at the source coordinates
# of feature X, so a drawn-coordinate hash of Y equals the source hash of X.
CROP = "201-4000"


def _linear(*inputs: RecordInput, **options) -> LinearDiagramRequest:
    options.setdefault("selected_features_set", TYPES)
    return LinearDiagramRequest(records=inputs, options=LinearDiagramOptions(**options))


def _collide(record_key: str = "k", region: str | None = CROP) -> RecordInput:
    return RecordInput(
        GenBankInputSource(COLLIDE),
        record_key=record_key,
        region=parse_region_spec(region) if region else None,
    )


def _source_id(plan, record_index: int, feature_type: str, start: int) -> str:
    """Biological feature ID of the source feature that starts at 1-based ``start``."""
    return next(
        entry.biological_feature_id
        for entry in plan.provenance[record_index].source_feature_catalog
        if entry.feature_type == feature_type
        and min(part[0] for part in entry.location_parts) + 1 == start
    )


def _drawn_id(plan, record_index: int, biological_feature_id: str) -> str:
    identity = FeatureIdentity(plan.provenance[record_index].record_key, biological_feature_id)
    binding = resolve_feature_identities(
        records=plan.records,
        record_keys=tuple(item.record_key for item in plan.provenance),
        source_catalogs=tuple(item.source_feature_catalog for item in plan.provenance),
        identities=(identity,),
    )[identity]
    return get_feature_hash(binding.feature, plan.records[record_index].id)


def _drawn(svg: str) -> set[tuple[str, int]]:
    return {
        (feature_id, int(index))
        for feature_id, index in re.findall(
            r'data-gbdraw-feature-id="([^"]+)"[^>]*data-gbdraw-record-index="(\d+)"', svg
        )
    }


def _labels(svg: str) -> list[tuple[str, str]]:
    return re.findall(r'<text[^>]*data-label-feature-id="([^"]+)"[^>]*>([^<]*)</text>', svg)


def _svg(request) -> str:
    return build_request_diagram(request).drawing.tostring()


def test_feature_off_hides_its_feature_and_not_the_one_drawn_at_its_source_coordinates():
    base = plan_request(_linear(_collide()))
    x, y = _source_id(base, 0, "misc_feature", 2401), _source_id(base, 0, "misc_feature", 2601)
    # The collision this identity avoids: Y is drawn where X was in the source.
    assert _drawn_id(base, 0, y) == x.split("~")[0]
    request = _linear(_collide(), feature_overrides=(FeatureOverride("k", x, feature_visibility="off"),))
    drawn = _drawn(_svg(request))
    assert (_drawn_id(base, 0, x), 0) not in drawn
    assert (_drawn_id(base, 0, y), 0) in drawn
    assert plan_request(request).inputs.feature_identity_notices == ()


def test_label_text_reaches_only_its_feature():
    base = plan_request(_linear(_collide()))
    tx = _source_id(base, 0, "tRNA", 2550)
    request = _linear(
        _collide(),
        feature_overrides=(FeatureOverride("k", tx, label_text="TX label"),),
        config_overrides={"labels.linear.scope": "all"},
    )
    labels = _labels(_svg(request))
    assert [text for _id, text in labels if text in {"TX label", "tRNA-Leu"}].count("TX label") == 1
    assert ("tRNA-Leu" in [text for _id, text in labels])
    assert (_drawn_id(base, 0, tx), "TX label") in labels


def test_duplicate_record_edits_apply_to_their_copy_only():
    request = CircularDiagramRequest(
        records=(RecordInput(
            GenBankInputSource(FIXTURES / "b_dup_ids.gb"), record_key="dup",
            cardinality=RecordCardinality.ALL,
        ),),
        options=CircularDiagramOptions(),
        layout=CircularMultiRecordOptions(),
        grouping="grid",
    )
    base = plan_request(request)
    first = _source_id(base, 0, "CDS", 301)
    edited = CircularDiagramRequest(
        records=request.records,
        options=CircularDiagramOptions(
            feature_overrides=(FeatureOverride("dup:1", first, feature_visibility="off"),)
        ),
        layout=request.layout,
        grouping="grid",
    )
    drawn = _drawn(_svg(edited))
    drawn_id = _drawn_id(base, 0, first)
    assert (drawn_id, 0) not in drawn
    assert (drawn_id, 1) in drawn


def _identical_features_record() -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 40), id="dup", annotations={"molecule_type": "DNA", "topology": "linear"}
    )
    record.features = [
        SeqFeature(SimpleLocation(20, 40, 1), type="CDS", qualifiers={"locus_tag": [tag]})
        for tag in ("first", "second")
    ]
    return record


def test_label_on_shows_only_the_named_instance_of_identical_features():
    request = _linear(RecordInput(InMemoryRecordSource(_identical_features_record()), record_key="one"))
    first, second = (entry.biological_feature_id for entry in plan_request(request).provenance[0].source_feature_catalog)
    assert second.endswith("~1")
    edited = _linear(
        RecordInput(InMemoryRecordSource(_identical_features_record()), record_key="one"),
        feature_overrides=(FeatureOverride("one", second, label_visibility="on"),),
        config_overrides={"labels.linear.scope": "none"},
    )
    assert [text for _id, text in _labels(_svg(edited))] == ["second"]


@pytest.mark.parametrize("scope,expected", [("none", []), ("all", ["Renamed"])])
def test_text_only_edit_changes_a_shown_label_and_never_shows_one(scope, expected):
    base = plan_request(_linear(_collide(region=None)))
    cds = _source_id(base, 0, "CDS", 2001)
    request = _linear(
        _collide(region=None),
        selected_features_set=("CDS",),
        feature_overrides=(FeatureOverride("k", cds, label_text="Renamed"),),
        config_overrides={"labels.linear.scope": scope},
    )
    texts = [text for _id, text in _labels(_svg(request))]
    assert [text for text in texts if text in {"Renamed", "gtg start"}] == expected


@pytest.mark.parametrize("scope,other_labels", [("none", 0), ("outer", 1)])
def test_circular_label_on_shows_its_label_under_any_scope(scope, other_labels):
    source = RecordInput(GenBankInputSource(COLLIDE), record_key="k")
    options = {"selected_features_set": ("CDS", "tRNA")}
    base = plan_request(CircularDiagramRequest(records=(source,), options=CircularDiagramOptions(**options)))
    tx = _source_id(base, 0, "tRNA", 2550)
    svg = _svg(CircularDiagramRequest(records=(source,), options=CircularDiagramOptions(
        **options,
        feature_overrides=(FeatureOverride("k", tx, label_visibility="on", label_text="TX label"),),
        config_overrides={"labels.circular.scope": scope},
    )))
    assert (svg.count(">TX label<"), svg.count("tRNA-Leu")) == (1, other_labels)


def test_label_on_without_text_falls_back_to_type_and_source_selector_location():
    # The same text as the Web's resolveDefaultLabelText: "<type> <selector location>".
    base = plan_request(_linear(_collide()))
    x = _source_id(base, 0, "misc_feature", 2401)
    request = _linear(
        _collide(),
        feature_overrides=(FeatureOverride("k", x, label_visibility="on"),),
        config_overrides={"labels.linear.scope": "none"},
    )
    assert [text for _id, text in _labels(_svg(request))] == ["misc_feature 2400..2500"]


def test_identity_annotation_attaches_to_its_feature_after_crop():
    base = plan_request(_linear(_collide()))
    x = _source_id(base, 0, "misc_feature", 2401)
    annotations = AnnotationOptions(sets=(AnnotationSet(
        id="marks",
        annotations=(RegionAnnotation(id="x", target=FeatureIdentitySpan("k", x), label="X"),),
    ),))
    plan = plan_request(_linear(_collide(), annotations=annotations))
    assert [item.segments for item in plan.resolved_annotations.annotations] == [((2200, 2300),)]
    assert plan.resolved_annotations.warnings == ()


def test_identity_annotation_outside_the_crop_is_skipped_as_unmatched():
    base = plan_request(_linear(_collide(region=None)))
    early = _source_id(base, 0, "CDS", 301)
    annotations = AnnotationOptions(sets=(AnnotationSet(
        id="marks", annotations=(RegionAnnotation(id="a", target=FeatureIdentitySpan("k", early)),),
    ),))
    plan = plan_request(_linear(_collide(region="2001-4000"), annotations=annotations))
    assert plan.resolved_annotations.annotations == ()
    assert [(w.code, w.missing_count) for w in plan.resolved_annotations.warnings] == [
        ("feature_selector_unmatched", 1)
    ]


def test_edits_that_are_not_drawn_notify_and_generate_succeeds():
    base = plan_request(_linear(_collide(region=None)))
    early = _source_id(base, 0, "CDS", 301)
    request = _linear(
        _collide(region="2001-4000"),
        feature_overrides=(
            FeatureOverride("k", early, feature_visibility="on", label_text="Early"),
            FeatureOverride("k", "fdeadbeef", label_visibility="off"),
        ),
    )
    plan = plan_request(request)
    assert [
        (n.biological_feature_id, n.status, n.kinds) for n in plan.inputs.feature_identity_notices
    ] == [
        (early, "crop_excluded", ("feature_visibility", "label_text")),
        ("fdeadbeef", "unresolved", ("label_visibility",)),
    ]
    assert build_request_diagram(request).feature_identity_notices == plan.inputs.feature_identity_notices


def test_identity_rows_decide_before_visibility_table_rules():
    base = plan_request(_linear(_collide()))
    x, y = _source_id(base, 0, "misc_feature", 2401), _source_id(base, 0, "misc_feature", 2601)
    table = DataFrame(
        [["*", "misc_feature", "location", ".*", "off"]],
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    )
    request = _linear(
        _collide(),
        feature_visibility_table=table,
        feature_overrides=(FeatureOverride("k", x, feature_visibility="on"),),
    )
    drawn = _drawn(_svg(request))
    assert (_drawn_id(base, 0, x), 0) in drawn
    assert (_drawn_id(base, 0, y), 0) not in drawn


def test_legend_lists_the_type_of_a_feature_shown_by_its_edit():
    base = plan_request(_linear(_collide()))
    x = _source_id(base, 0, "misc_feature", 2401)

    def legend(**options):
        svg = _svg(_linear(_collide(), selected_features_set=("CDS",), **options))
        return {text for text in re.findall(r">([^<>]+)</text>", svg) if text in {"CDS", "misc_feature"}}

    assert legend() == {"CDS"}
    assert legend(feature_overrides=(FeatureOverride("k", x, feature_visibility="on"),)) == {
        "CDS", "misc_feature"
    }


def test_exclude_matching_keeps_the_type_selection_and_skips_table_rules():
    base = plan_request(_linear(_collide()))
    x, y = _source_id(base, 0, "misc_feature", 2401), _source_id(base, 0, "misc_feature", 2601)
    table = DataFrame(
        [["*", "misc_feature", "location", ".*", "show"]],
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    )
    request = _linear(
        _collide(),
        selected_features_set=("CDS",),
        feature_visibility_table=table,
        feature_overrides=(FeatureOverride("k", x, feature_visibility="exclude_matching"),),
    )
    drawn = _drawn(_svg(request))
    assert (_drawn_id(base, 0, x), 0) not in drawn
    assert (_drawn_id(base, 0, y), 0) in drawn


def _protein_tags(extraction) -> list[str]:
    return [protein.locus_tag for protein in extraction.proteins_by_record[0]]


@pytest.mark.parametrize("mode", ["off", "exclude_matching"])
def test_hidden_and_excluded_cds_leave_losatp_protein_extraction(mode):
    base = plan_request(_linear(_collide(region=None)))
    cds = _source_id(base, 0, "CDS", 2001)
    request = _linear(
        _collide(region=None),
        feature_overrides=(FeatureOverride("k", cds, feature_visibility=mode),),
    )
    plan = plan_request(request)
    before = _protein_tags(_extract_linear_request_proteins(base.request, base.records, base.inputs))
    after = _protein_tags(_extract_linear_request_proteins(plan.request, plan.records, plan.inputs))
    assert "TESTA_0006" in before
    assert after == [tag for tag in before if tag != "TESTA_0006"]


def test_browser_protein_helper_resolves_rows_with_the_shared_resolver():
    from tests.test_protein_colinearity import _load_web_helper_namespace

    namespace = _load_web_helper_namespace()
    base = plan_request(_linear(_collide(region=None)))
    cds = _source_id(base, 0, "CDS", 2001)

    def extract(rows):
        result = json.loads(str(namespace["extract_cds_protein_fasta"](
            str(COLLIDE), "genbank", None, None, None, "1", 0, "k", None,
            None if rows is None else json.dumps(rows),
        )))
        assert "error" not in result, result
        return sorted(item["locus_tag"] for item in result["protein_map"].values())

    row = FeatureOverride("k", cds, feature_visibility="off").to_mapping()
    assert "TESTA_0006" in extract(None)
    assert extract([row]) == [tag for tag in extract(None) if tag != "TESTA_0006"]


def test_gff_feature_shown_by_its_edit_is_loaded_and_drawn(tmp_path):
    gff, fasta = tmp_path / "source.gff3", tmp_path / "source.fasta"
    fasta.write_text(">record\n" + "ATG" * 40 + "\n")
    gff.write_text(
        "##gff-version 3\n"
        "record\t.\tgene\t1\t18\t.\t+\t.\tID=hidden-gene\n"
        "record\t.\tCDS\t31\t60\t.\t+\t0\tID=shown-cds\n"
    )
    source = RecordInput(GffFastaInputSource(gff, fasta), record_key="one")
    gene = plan_request(
        _linear(source, selected_features_set=("CDS",))
    ).provenance[0].source_feature_catalog[0].biological_feature_id
    request = _linear(
        source,
        selected_features_set=("CDS",),
        feature_overrides=(FeatureOverride("one", gene, feature_visibility="on"),),
    )
    plan = plan_request(request)
    assert plan.inputs.feature_identity_notices == ()
    assert plan.inputs.record_features[0].overrides
    assert {feature.type for feature in plan.records[0].features} >= {"gene", "CDS"}


def test_package_api_reports_notices_for_unresolved_placements():
    # Owner Q3 = A: an exact placement whose identity the source lacks notifies.
    records = read_genbank(COLLIDE)
    diagram = draw_linear(
        records,
        options=LinearOptions(features=FeatureOptions(
            types=TYPES,
            placements=(FeaturePlacementOverride("record-1", "fdeadbeef", FeaturePlacementTarget("main")),),
        )),
    )
    assert [(n.record_key, n.status) for n in diagram.feature_identity_notices] == [
        ("record-1", "unresolved")
    ]


def test_cli_session_replay_applies_rows_and_logs_notices(tmp_path, caplog):
    from gbdraw.linear import linear_main
    from gbdraw.session import save_session_document

    # Python session writers store the resolved (cropped) record, so this uses the whole record.
    base = plan_request(_linear(_collide(region=None)))
    x, y = _source_id(base, 0, "misc_feature", 2401), _source_id(base, 0, "misc_feature", 2601)
    request = _linear(
        _collide(region=None),
        feature_overrides=(
            FeatureOverride("k", x, feature_visibility="off"),
            FeatureOverride("k", "fdeadbeef", feature_visibility="off"),
        ),
    )
    session = tmp_path / "identity.gbdraw-session.json"
    document = save_session_document(session, request)
    assert document.to_dict()["renderRequest"]["diagramOptions"]["featureOverrides"] == [
        row.to_mapping() for row in request.options.feature_overrides
    ]
    prefix = tmp_path / "replayed"
    with caplog.at_level(logging.WARNING):
        linear_main(["--session", str(session), "-o", str(prefix), "-f", "svg"])
    drawn = _drawn(prefix.with_suffix(".svg").read_text(encoding="utf-8"))
    assert (_drawn_id(base, 0, x), 0) not in drawn
    assert (_drawn_id(base, 0, y), 0) in drawn
    assert any(
        "feature_identity_unresolved" in message and "fdeadbeef" in message
        for message in caplog.messages
    )


def test_override_rows_reject_empty_or_multiline_edits():
    from gbdraw.exceptions import ValidationError

    with pytest.raises(ValidationError, match="at least one edit"):
        FeatureOverride("k", "f1")
    with pytest.raises(ValidationError, match="one non-blank line"):
        FeatureOverride("k", "f1", label_text="two\nlines")
    with pytest.raises(ValidationError, match="Duplicate feature override identity"):
        LinearDiagramOptions(feature_overrides=(
            FeatureOverride("k", "f1", label_text="a"), FeatureOverride("k", "f1", label_visibility="on"),
        ))


def test_collide_fixture_layout_is_as_documented():
    record = next(SeqIO.parse(COLLIDE, "genbank"))
    spans = [(f.type, int(f.location.start) + 1, int(f.location.end)) for f in record.features if f.type in {"misc_feature", "tRNA"}]
    assert spans == [
        ("misc_feature", 2401, 2500), ("tRNA", 2550, 2620),
        ("misc_feature", 2601, 2700), ("tRNA", 2750, 2820),
    ]
