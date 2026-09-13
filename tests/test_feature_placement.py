"""Requested placement identity; geometry belongs to the following planner phase."""

from dataclasses import FrozenInstanceError

import pytest


from dataclasses import asdict, replace
import copy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.prepared import PreparedBiologicalInputCache, PreparedResourceIdentity
from gbdraw.api.request_render import plan_request, resolve_request
from gbdraw.api.requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    GenBankInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordDisplayOptions,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.placement import FeaturePlacementOverride, FeaturePlacementTarget
from gbdraw.features.source import build_source_feature_catalog
from gbdraw.io.regions import parse_region_spec
from gbdraw.session_request_codec import (
    CanonicalRequestEncodingError,
    encode_canonical_request,
)


def test_main_intent_is_immutable_and_has_no_transient_lane():
    from gbdraw.features.placement import (
        FeaturePlacementOverride,
        FeaturePlacementTarget,
    )

    target = FeaturePlacementTarget(kind="main")
    override = FeaturePlacementOverride("record-1", "f1234abcd", target)
    assert override.target.kind == "main"
    assert override.target.side is None and override.target.level is None
    assert not hasattr(override, "feature_track_id")
    with pytest.raises(FrozenInstanceError):
        override.record_key = "record-2"


def _record():
    record = SeqRecord(
        Seq("ATG" * 40),
        id="duplicate",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
        },
    )
    record.features = [
        SeqFeature(
            SimpleLocation(2, 8, 1), type="gene", qualifiers={"locus_tag": ["early"]}
        ),
        SeqFeature(
            SimpleLocation(20, 40, 1), type="CDS", qualifiers={"locus_tag": ["alpha"]}
        ),
        SeqFeature(
            SimpleLocation(20, 40, 1), type="CDS", qualifiers={"locus_tag": ["beta"]}
        ),
        SeqFeature(
            CompoundLocation([SimpleLocation(60, 75, -1), SimpleLocation(45, 55, -1)]),
            type="CDS",
            qualifiers={"locus_tag": ["multipart"]},
        ),
        SeqFeature(
            SimpleLocation(85, 100, -1),
            type="repeat_region",
            qualifiers={"locus_tag": ["repeat"]},
        ),
    ]
    return record


def _request(record=None, *, mode="circular", records=None, **options):
    cls, opts = (
        (CircularDiagramRequest, CircularDiagramOptions)
        if mode == "circular"
        else (LinearDiagramRequest, LinearDiagramOptions)
    )
    return cls(
        records=tuple(
            records
            or (
                RecordInput(
                    InMemoryRecordSource(record or _record()), record_key="one"
                ),
            )
        ),
        options=opts(**options),
    )


def _table(selector="locus_tag=alpha", placement="main", **extra):
    return DataFrame([dict(feature_selector=selector, placement=placement, **extra)])


def _exact(record, index=1, key="one", target=None):
    return FeaturePlacementOverride(
        key,
        build_source_feature_catalog(record)[index].biological_feature_id,
        target or FeaturePlacementTarget("main"),
    )


@pytest.mark.parametrize(
    "mode,side",
    [
        ("circular", "outward"),
        ("circular", "inward"),
        ("linear", "above"),
        ("linear", "below"),
    ],
)
def test_directional_shape_and_deferred_geometry_validation(mode, side):
    target = FeaturePlacementTarget.from_mapping(
        {"kind": "lane", "side": side, "level": 1}
    )
    override = _exact(_record(), target=target)
    # Separate strands and one-sided slot presets require the later resolved-slot planner.
    config = {"canvas.strandedness": True, "canvas.resolve_overlaps": False}
    plan = plan_request(
        _request(mode=mode, feature_placements=(override,), config_overrides=config)
    )
    assert plan.inputs.placements[0].foreground[0].target == target
    assert asdict(target) == {"kind": "lane", "side": side, "level": 1}
    assert not hasattr(plan.inputs.placements[0].foreground[0], "feature_track_id")


@pytest.mark.parametrize(
    "payload",
    [
        None,
        [],
        {},
        {"kind": "main", "side": None},
        {"kind": "main", "level": 1},
        {"kind": "auto"},
        {"kind": "lane"},
        {"kind": "main", "pixel": 5},
        {"kind": "lane", "side": "outward", "level": 1, "feature_track_id": 1},
        {"kind": "lane", "side": "unknown", "level": 1},
        *[
            {"kind": "lane", "side": "outward", "level": v}
            for v in (0, 2, -1, True, 1.0, "1", None)
        ],
    ],
)
def test_invalid_target_shape_tokens_levels(payload):
    with pytest.raises(ValidationError):
        FeaturePlacementTarget.from_mapping(payload)


@pytest.mark.parametrize(
    "payload",
    [
        None,
        {},
        {
            "recordKey": "one",
            "biologicalFeatureId": "f1234abcd",
            "placement": {"kind": "main"},
            "svgId": "f1",
        },
    ],
)
def test_exact_unknown_fields_rejected(payload):
    with pytest.raises(ValidationError):
        FeaturePlacementOverride.from_mapping(payload)


@pytest.mark.parametrize(
    "record_key,feature_id",
    [("", "f1"), ("one\0f1", "f1"), ("one", "f1\0"), (None, "f1")],
)
def test_invalid_exact_identity(record_key, feature_id):
    with pytest.raises(ValidationError):
        FeaturePlacementOverride(record_key, feature_id, FeaturePlacementTarget("main"))


@pytest.mark.parametrize("mode,side", [("circular", "above"), ("linear", "inward")])
def test_wrong_mode_exact_and_table_rejected(mode, side):
    with pytest.raises(ValidationError, match="unsupported.*mode"):
        _request(
            mode=mode,
            feature_placements=(
                _exact(_record(), target=FeaturePlacementTarget("lane", side, 1)),
            ),
        )
    with pytest.raises(ValidationError, match="unsupported.*mode"):
        plan_request(
            _request(mode=mode, feature_placement_table=_table(placement=side))
        )


@pytest.mark.parametrize("kind", ["duplicate", "table", "file", "table_file"])
def test_duplicate_or_mixed_input_rejected(kind, tmp_path):
    exact = _exact(_record())
    kwargs = {"feature_placements": (exact,)}
    if kind == "duplicate":
        kwargs["feature_placements"] = (exact, exact)
    elif kind == "table":
        kwargs["feature_placement_table"] = _table()
    elif kind == "file":
        kwargs["feature_placement_table_file"] = tmp_path / "placements.tsv"
    else:
        kwargs = dict(
            feature_placement_table=_table(),
            feature_placement_table_file="placements.tsv",
        )
    with pytest.raises(ValidationError, match="Duplicate|mutually exclusive"):
        _request(**kwargs)


@pytest.mark.parametrize("selector,count", [("locus_tag=missing", 0), ("type=CDS", 3)])
def test_source_selector_zero_and_many_error(selector, count):
    with pytest.raises(ValidationError, match=f"matched {count}"):
        plan_request(_request(feature_placement_table=_table(selector)))


def test_selector_normalizes_auto_removes_table_and_keeps_original_dataframe():
    table = DataFrame(
        [
            dict(
                record="#1",
                feature_selector="locus_tag=beta",
                placement=" Outward ",
                level="",
            ),
            dict(
                record="#1",
                feature_selector="locus_tag=alpha",
                placement="AUTO",
                level="",
            ),
        ]
    )
    before = table.copy(deep=True)
    request = _request(feature_placement_table=table)
    plan = plan_request(request)
    exact = plan.request.options.feature_placements
    assert len(exact) == 1 and exact[0].target == FeaturePlacementTarget(
        "lane", "outward", 1
    )
    assert exact[0].biological_feature_id.endswith("~2")
    assert plan.request.options.feature_placement_table is None
    assert plan.request.options.feature_placement_table_file is None
    assert request.options.feature_placement_table is table and table.equals(before)
    assert plan.inputs.placements[0].overrides[0].source_feature_index == 2


def test_file_table_supports_hash_record_index_and_drops_path(tmp_path):
    path = tmp_path / "placement.tsv"
    path.write_text(
        "record\tfeature_selector\tplacement\tlevel\n#1\tlocus_tag=alpha\tmain\t\n",
        encoding="utf-8",
    )
    plan = plan_request(_request(feature_placement_table_file=path))
    assert plan.request.options.feature_placements == (_exact(_record()),)
    assert plan.request.options.feature_placement_table_file is None


@pytest.mark.parametrize(
    "rows",
    [
        [dict(feature_selector="locus_tag=alpha", placement="main", unknown="x")],
        [dict(feature_selector="locus_tag=alpha", placement="main", level="1")],
        [dict(feature_selector="locus_tag=alpha", placement="auto", level="1")],
        [dict(feature_selector="locus_tag=alpha", placement="below", level=True)],
        [dict(feature_selector="locus_tag=alpha", placement="outward", level="2")],
        [dict(feature_selector="locus_tag=alpha", placement="pixel")],
        [
            dict(feature_selector="locus_tag=alpha", placement=value)
            for value in ("main", "auto")
        ],
    ],
)
def test_invalid_table_rows(rows):
    with pytest.raises(ValidationError):
        plan_request(_request(feature_placement_table=DataFrame(rows)))


def test_dataframe_optional_levels_preserve_integer_lane_one():
    # A numeric level mixed with blank Main/Auto cells becomes float64 in pandas.
    table = DataFrame([
        dict(feature_selector="locus_tag=alpha", placement="above", level=1),
        dict(feature_selector="locus_tag=early", placement="main", level=None),
        dict(feature_selector="locus_tag=beta", placement="auto", level=None),
    ])
    plan = plan_request(_request(mode="linear", feature_placement_table=table))
    targets = {item.target for item in plan.request.options.feature_placements}
    assert targets == {FeaturePlacementTarget("main"), FeaturePlacementTarget("lane", "above", 1)}
    assert plan.request.options.feature_placement_table is None


def test_exact_sort_and_duplicate_hashes_are_independent():
    record = _record()
    alpha, beta = _exact(record), _exact(record, 2)
    assert alpha.biological_feature_id != beta.biological_feature_id
    plan = plan_request(_request(record, feature_placements=(beta, alpha)))
    assert plan.request.options.feature_placements == (alpha, beta)
    assert [
        row.source_feature_index for row in plan.inputs.placements[0].foreground
    ] == [1, 2]


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("start", [None, 1, 25])
def test_source_identity_intent_and_parts_survive_display_and_reverse(
    mode, reverse, start
):
    record = _record()
    before = copy.deepcopy(record)
    exact = _exact(record, 3)
    inp = RecordInput(
        InMemoryRecordSource(record),
        record_key="one",
        display=RecordDisplayOptions(start_coordinate=start),
        presentation=RecordPresentation(reverse_complement=reverse),
    )
    plan = plan_request(
        _request(mode=mode, records=(inp,), feature_placements=(exact,))
    )
    assert plan.request.options.feature_placements == (exact,)
    binding = plan.inputs.placements[0].foreground[0]
    assert (
        binding.source_feature_index == 3
        and binding.biological_feature_id == exact.biological_feature_id
    )
    assert plan.provenance[0].source_feature_catalog[3].location_parts == (
        (60, 75, -1),
        (45, 55, -1),
    )
    assert (
        str(record.seq) == str(before.seq) and record.annotations == before.annotations
    )
    assert [vars(f) for f in record.features] == [vars(f) for f in before.features]


@pytest.mark.parametrize(
    "region,status", [("25-70", "foreground"), ("80-120", "crop_excluded")]
)
@pytest.mark.parametrize("reverse", [False, True])
def test_request_crop_known_dormant_and_materialized_request_rebind(
    region, status, reverse
):
    record = _record()
    exact = _exact(record)
    if reverse:
        start, end = map(int, region.split("-"))
        region = f"{len(record) + 1 - end}-{len(record) + 1 - start}"
    request = _request(
        records=(
            RecordInput(
                InMemoryRecordSource(record),
                record_key="one",
                region=parse_region_spec(region),
                presentation=RecordPresentation(reverse_complement=reverse),
            ),
        ),
        feature_placements=(exact,),
    )
    plan = plan_request(request)
    assert plan.inputs.placements[0].overrides[0].status == status
    assert bool(plan.inputs.placements[0].foreground) == (status == "foreground")
    replay = plan_request(plan.request)
    assert replay.inputs.placements == plan.inputs.placements
    assert (
        replay.provenance[0].source_feature_catalog
        is plan.provenance[0].source_feature_catalog
    )
    restored = plan_request(
        replace(request, records=(replace(request.records[0], region=None),))
    )
    assert (
        restored.inputs.placements[0].foreground[0].biological_feature_id
        == exact.biological_feature_id
    )


def test_selector_uses_original_source_before_crop():
    request = _request(
        records=(
            RecordInput(
                InMemoryRecordSource(_record()),
                record_key="one",
                region=parse_region_spec("80-120"),
            ),
        ),
        feature_placement_table=_table(),
    )
    plan = plan_request(request)
    assert plan.inputs.placements[0].overrides[0].status == "crop_excluded"
    assert len(plan.request.options.feature_placements) == 1


@pytest.mark.parametrize("kind", ["hidden", "underlay", "foreground"])
def test_existing_visibility_and_rendering_owners_classify_occupancy(kind):
    record = _record()
    exact = _exact(record)
    options = {"feature_placements": (exact,)}
    if kind == "hidden":
        options["feature_visibility_table"] = DataFrame(
            [
                dict(
                    record_id="*",
                    feature_type="CDS",
                    qualifier="locus_tag",
                    value="^alpha$",
                    action="hide",
                )
            ]
        )
    elif kind == "underlay":
        options["feature_shapes"] = {"CDS": "underlay"}
    plan = plan_request(_request(record, **options))
    assert plan.inputs.placements[0].overrides[0].status == kind
    assert bool(plan.inputs.placements[0].foreground) == (kind == "foreground")
    assert plan.request.options.feature_placements == (exact,)
    restored = plan_request(
        replace(
            plan.request,
            options=replace(
                plan.request.options, feature_shapes=None, feature_visibility_table=None
            ),
        )
    )
    assert (
        restored.inputs.placements[0].foreground[0].biological_feature_id
        == exact.biological_feature_id
    )


def test_default_underlay_and_type_hidden_are_dormant():
    record = _record()
    request = _request(
        record,
        selected_features_set=("CDS", "repeat_region"),
        feature_placements=(_exact(record, 0), _exact(record, 4)),
    )
    statuses = {
        row.source_feature_index: row.status
        for row in plan_request(request).inputs.placements[0].overrides
    }
    assert statuses == {0: "hidden", 4: "underlay"}


@pytest.mark.parametrize(
    "key,identity", [("one", "stale"), ("missing", "stale"), ("one", "f1234_record_1")]
)
def test_unknown_or_rendered_identity_never_becomes_dormant(key, identity):
    override = FeaturePlacementOverride(key, identity, FeaturePlacementTarget("main"))
    with pytest.raises(ValidationError, match="Unknown|stale"):
        plan_request(_request(feature_placements=(override,)))


def test_caller_precropped_record_does_not_infer_unseen_source():
    original = _record()
    exact = _exact(original)
    crop = plan_request(
        _request(
            records=(
                RecordInput(
                    InMemoryRecordSource(original),
                    record_key="one",
                    region=parse_region_spec("80-120"),
                ),
            )
        )
    ).records[0]
    with pytest.raises(ValidationError, match="Unknown/stale"):
        plan_request(_request(crop, feature_placements=(exact,)))


@pytest.mark.parametrize("reorder", [False, True])
def test_duplicate_biological_record_ids_are_isolated_and_reordered_by_key(reorder):
    record = _record()
    records = [
        RecordInput(InMemoryRecordSource(record), record_key=key)
        for key in ("first", "second")
    ]
    if reorder:
        records.reverse()
    exact = _exact(record, key="second")
    plan = plan_request(
        _request(mode="linear", records=records, feature_placements=(exact,))
    )
    assert [r.record_key for r in plan.inputs.placements] == [
        r.record_key for r in records
    ]
    assert {r.record_key: len(r.foreground) for r in plan.inputs.placements} == {
        "first": 0,
        "second": 1,
    }
    assert (
        plan.provenance[0].source_feature_catalog
        is plan.provenance[1].source_feature_catalog
    )
    with pytest.raises(ValidationError, match="multiple records"):
        plan_request(
            _request(
                mode="linear",
                records=records,
                feature_placement_table=_table(record="duplicate"),
            )
        )
    selected = plan_request(
        _request(
            mode="linear", records=records, feature_placement_table=_table(record="#2")
        )
    )
    assert (
        selected.request.options.feature_placements[0].record_key
        == records[1].record_key
    )


def test_resource_catalog_is_built_once_and_reused_across_crop_reverse_display(
    tmp_path, monkeypatch
):
    import gbdraw.api.record_planning as owner

    path = tmp_path / "source.gbk"
    SeqIO.write(_record(), path, "genbank")
    cache = PreparedBiologicalInputCache()
    identity = PreparedResourceIdentity("source", "stable", path.stat().st_size)
    calls = []
    build = owner.build_source_feature_catalog
    monkeypatch.setattr(
        owner,
        "build_source_feature_catalog",
        lambda r: (calls.append(r.id), build(r))[1],
    )
    catalog = None
    for start, crop, reverse in [
        (None, None, False),
        (25, None, False),
        (None, "25-70", True),
    ]:
        inp = RecordInput(
            GenBankInputSource(path),
            record_key="one",
            region=parse_region_spec(crop) if crop else None,
            display=RecordDisplayOptions(start_coordinate=start),
            presentation=RecordPresentation(reverse_complement=reverse),
        )
        diagnostics = {"metrics": {}}
        with cache.transaction(
            resource_paths={path: identity}, diagnostics=diagnostics
        ):
            plan = plan_request(
                _request(records=(inp,), feature_placement_table=_table())
            )
        if catalog is not None:
            assert plan.provenance[0].source_feature_catalog is catalog
        catalog = plan.provenance[0].source_feature_catalog
        layers = _planned_layers(plan)
        assert len(layers.foreground_features) == (1 if crop else 3)
        if crop:
            assert plan.inputs.placements[0].overrides[0].status == "crop_excluded"
            assert all(f.placement.requested_target is None for f in layers.foreground_features.values())
        else:
            assert _assignments(layers)[1].requested_target == FeaturePlacementTarget("main")
    assert calls == ["duplicate"]
    assert diagnostics["metrics"]["parsedSourceParseCount"] == 0


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_current_writer_encodes_exact_and_rejects_unmaterialized_placement(mode):
    payload = encode_canonical_request(_request(mode=mode, feature_placements=(_exact(_record()),))).payload
    assert len(payload["diagramOptions"]["featurePlacements"]) == 1
    with pytest.raises(CanonicalRequestEncodingError, match="placement tables must be materialized"):
        encode_canonical_request(_request(mode=mode, feature_placement_table=_table()))


def test_auto_only_materializes_empty_and_retains_current_writer():
    request = resolve_request(
        _request(feature_placement_table=_table(placement="auto"))
    )
    assert (
        request.options.feature_placements == ()
        and request.options.feature_placement_table is None
    )
    payload = encode_canonical_request(request)
    assert payload.payload["diagramOptions"]["featurePlacements"] == []


def test_source_catalog_is_immutable_and_only_requested_types_are_public():
    import gbdraw
    import gbdraw.api

    catalog = build_source_feature_catalog(_record())
    with pytest.raises(FrozenInstanceError):
        catalog[0].source_feature_index = 100
    assert isinstance(catalog[0].qualifiers, tuple)
    for namespace in (gbdraw, gbdraw.api):
        assert "FeaturePlacementTarget" in namespace.__all__
        assert "FeaturePlacementOverride" in namespace.__all__


def test_batch_projects_only_each_instances_intent_and_binding():
    record = _record()
    request = CircularBatchRequest(
        records=tuple(
            RecordInput(InMemoryRecordSource(record), record_key=key)
            for key in ("first", "second")
        ),
        outputs=tuple(
            RenderOutputRequest(output_prefix=key) for key in ("first", "second")
        ),
        options=CircularDiagramOptions(
            feature_placements=(_exact(record, key="second"),)
        ),
    )
    plan = plan_request(request)
    first, second = plan.item_plans()
    assert first.request.options.feature_placements == ()
    assert first.inputs.placements[0].record_key == "first"
    assert first.inputs.placements[0].overrides == ()
    assert second.request.options.feature_placements == (_exact(record, key="second"),)
    assert len(second.inputs.placements) == 1
    assert second.inputs.placements[0] == plan.inputs.placements[1]
    assert first.inputs.features is second.inputs.features is plan.inputs.features
    assert plan_request(second.request).inputs.placements == second.inputs.placements


def test_annotation_single_selector_still_allows_multiple_source_matches():
    from gbdraw.annotations import (
        AnnotationSet,
        FeatureSelector,
        FeatureSpan,
        RegionAnnotation,
        resolve_annotations,
    )

    record = _record()
    target = FeatureSpan(
        None, (FeatureSelector("CDS", "type"),), envelope="outer_bounds"
    )
    result = resolve_annotations(
        (AnnotationSet("s", (RegionAnnotation("a", target),)),), [record], mode="linear"
    )
    assert result.annotations[0].segments == ((20, 75),)
    with pytest.raises(ValidationError, match="exactly one"):
        plan_request(_request(record, feature_placement_table=_table("type=CDS")))


def test_historical_reader_does_not_admit_placement_fields(tmp_path):
    from gbdraw.session_request_codec import (
        CanonicalRequestDecodingError,
        decode_canonical_request,
    )

    encoded = encode_canonical_request(_request())
    for resource in encoded.resources:
        (tmp_path / resource.resource_id).write_bytes(resource.content)
    for field in (
        "featurePlacements",
        "featurePlacementTable",
        "featurePlacementTableFile",
    ):
        payload = copy.deepcopy(encoded.payload)
        payload["schema"] = 6
        payload["records"][0].pop("display")
        payload["diagramOptions"].pop("featurePlacements")
        payload["diagramOptions"][field] = []
        with pytest.raises(CanonicalRequestDecodingError, match="Unknown field"):
            decode_canonical_request(
                payload,
                resource_paths={
                    resource.resource_id: tmp_path / resource.resource_id
                    for resource in encoded.resources
                },
                output_directory=tmp_path,
            )


def test_gff_source_catalog_retains_features_before_loader_visibility_filter(tmp_path):
    from gbdraw.api.requests import GffFastaInputSource

    gff = tmp_path / "source.gff3"
    fasta = tmp_path / "source.fasta"
    fasta.write_text(">record\n" + "ATG" * 40 + "\n")
    gff.write_text(
        "##gff-version 3\n"
        "record\t.\tgene\t1\t18\t.\t+\t.\tID=hidden-gene\n"
        "record\t.\tCDS\t31\t60\t.\t+\t0\tID=shown-cds\n"
    )
    request = _request(
        records=(RecordInput(GffFastaInputSource(gff, fasta), record_key="one"),),
        selected_features_set=("CDS",),
        feature_placement_table=_table("ID=hidden-gene"),
    )
    plan = plan_request(request)
    assert len(plan.provenance[0].source_feature_catalog) == 2
    assert [feature.type for feature in plan.records[0].features] == ["CDS"]
    assert plan.inputs.placements[0].overrides[0].status == "hidden"
    assert plan.inputs.placements[0].foreground == ()
    shown = plan_request(
        replace(
            request,
            options=replace(request.options, selected_features_set=("CDS", "gene")),
        )
    )
    assert (
        shown.request.options.feature_placements
        == plan.request.options.feature_placements
    )
    assert shown.inputs.placements[0].foreground[0].source_feature_index == 0


def test_gff_catalog_capture_keeps_duplicate_biological_records_instance_aligned(
    tmp_path, monkeypatch
):
    import gbdraw.io.genome as genome
    from gbdraw.api.requests import GffFastaInputSource, RecordCardinality

    first, second = _record(), _record()
    first.features = first.features[:1]
    second.features = second.features[3:4]
    for record in (first, second):
        for feature in record.features:
            feature.sub_features = []
    monkeypatch.setattr(genome.GFF, "parse", lambda path: iter((first, second)))
    gff, fasta = tmp_path / "source.gff3", tmp_path / "source.fasta"
    gff.write_text("##gff-version 3\n")
    SeqIO.write(_record(), fasta, "fasta")
    request = _request(
        mode="linear",
        records=(
            RecordInput(
                GffFastaInputSource(gff, fasta),
                record_key="source",
                cardinality=RecordCardinality.ALL,
            ),
        ),
        feature_placement_table=_table("locus_tag=early", record="#1"),
    )
    plan = plan_request(request)
    assert [item.record_key for item in plan.provenance] == ["source:1", "source:2"]
    assert plan.provenance[0].source_feature_catalog[0].feature_type == "gene"
    assert plan.provenance[1].source_feature_catalog[0].feature_type == "CDS"
    assert plan.request.options.feature_placements[0].record_key == "source:1"


def test_shared_planner_reserves_fixed_main_before_auto():
    from gbdraw.features.factory import create_feature_layers
    from gbdraw.features.placement import FeaturePlacementSlot
    plan = plan_request(_request(feature_placement_table=_table("locus_tag=beta")))
    layers = create_feature_layers(
        plan.records[0], {}, ["CDS"], {"default": "#999999", "CDS": "#999999"},
        False, True, {}, compute_label_text=False,
        placement_inputs=plan.inputs.placements[0],
        placement_slot=FeaturePlacementSlot("circular", "split", False),
    )
    features = {f.source_feature_index: f for f in layers.foreground_features.values()}
    assert features[2].placement.level == 0
    assert features[2].placement.requested_target == FeaturePlacementTarget("main")
    assert (features[1].placement.side, features[1].placement.level) == ("outward", 1)
    assert features[3].placement.level == 0
    with pytest.raises(FrozenInstanceError):
        features[1].placement.level = 4


def _planned_layers(plan, *, mode="circular", direction="split", separate=False,
                    resolve=True, tolerance=0, index=0):
    from gbdraw.features.factory import create_feature_layers
    from gbdraw.features.placement import FeaturePlacementSlot
    return create_feature_layers(
        plan.records[index], {}, ["CDS", "repeat_region"],
        {"default": "#999999", "CDS": "#999999", "repeat_region": "#999999"},
        separate, resolve, {}, compute_label_text=False,
        feature_shapes=plan.request.options.feature_shapes,
        record_transform=plan.transforms[index],
        placement_inputs=plan.inputs.placements[index] if plan.inputs.placements else None,
        placement_slot=FeaturePlacementSlot(mode, direction, separate),
        feature_overlap_tolerance_bp=tolerance,
    )


def _assignments(layers):
    return {f.source_feature_index: f.placement for f in layers.foreground_features.values()}


@pytest.mark.parametrize("mode,direction", [("circular", "split"), ("linear", "overlay")])
@pytest.mark.parametrize("resolve", [False, True])
@pytest.mark.parametrize("fixed", [(), (1,), (2,), (0, 1, 2)])
@pytest.mark.parametrize("tolerance", [19, 20])
@pytest.mark.parametrize("unknown_strand", [None, 0])
def test_undefined_shares_negative_pool_for_fixed_and_auto(
    mode, direction, resolve, fixed, tolerance, unknown_strand,
):
    # Approved bounded legacy change: undefined and negative occupy one physical
    # pool. Positive Main remains independent even at the same coordinates.
    source = _record()
    source.features = [
        SeqFeature(SimpleLocation(10, 90, 1), type="CDS", qualifiers={"locus_tag": ["positive"]}),
        SeqFeature(SimpleLocation(20, 50, -1), type="CDS", qualifiers={"locus_tag": ["negative"]}),
        SeqFeature(SimpleLocation(30, 50, unknown_strand), type="CDS",
                   qualifiers={"locus_tag": ["undefined"]}),
    ]
    plan = plan_request(_request(
        source, mode=mode, feature_placements=tuple(_exact(source, i) for i in fixed),
    ))
    if 1 in fixed and 2 in fixed and tolerance < 20:
        with pytest.raises(ValidationError, match="Fixed feature placement conflict"):
            _planned_layers(plan, mode=mode, direction=direction, separate=True,
                            resolve=resolve, tolerance=tolerance)
        return
    layers = _planned_layers(plan, mode=mode, direction=direction, separate=True,
                             resolve=resolve, tolerance=tolerance)
    values = _assignments(layers)
    levels = [0, 0, 0]
    if resolve and tolerance < 20:
        levels[1 if 2 in fixed else 2] = 1
    assert [(values[i].strand_pool, values[i].level) for i in range(3)] == [
        ("positive", levels[0]), ("negative", levels[1]), ("negative", levels[2]),
    ]
    assert [values[i].requested_target for i in range(3)] == [
        FeaturePlacementTarget("main") if i in fixed else None for i in range(3)
    ]
    assert {values[i].side for i in (1, 2)} == {"inward" if mode == "circular" else "below"}
    assert {f.source_feature_index: f.feature_track_id for f in layers.foreground_features.values()} == {
        0: 0, 1: -1-levels[1], 2: -1-levels[2],
    }


@pytest.mark.parametrize("mode,direction", [
    ("circular", "split"), ("circular", "inside"), ("circular", "outside"),
    ("linear", "overlay"), ("linear", "above"), ("linear", "below"),
])
@pytest.mark.parametrize("separate", [False, True])
@pytest.mark.parametrize("resolve", [False, True])
def test_resolved_slot_main_and_auto_matrix(mode, direction, separate, resolve):
    plan = plan_request(_request(mode=mode, feature_placement_table=_table("locus_tag=beta")))
    result = _assignments(_planned_layers(
        plan, mode=mode, direction=direction, separate=separate, resolve=resolve,
    ))
    assert result[2].level == 0
    assert result[2].requested_target == FeaturePlacementTarget("main")
    assert result[1].level == int(resolve)
    assert result[3].level == 0
    assert result[3].strand_pool == ("negative" if separate else "combined")
    expected_side = {
        ("circular", "inside"): "inward", ("circular", "outside"): "outward",
        ("linear", "above"): "above", ("linear", "below"): "below",
    }.get((mode, direction), ("inward" if mode == "circular" else "below") if separate else "main")
    assert result[3].side == expected_side


@pytest.mark.parametrize("mode,direction,side", [
    ("circular", "split", "outward"), ("circular", "split", "inward"),
    ("linear", "overlay", "above"), ("linear", "overlay", "below"),
])
@pytest.mark.parametrize("resolve", [False, True])
def test_fixed_directional_resolver_matrix_and_distinct_physical_lanes(mode, direction, side, resolve):
    rows = DataFrame([
        dict(feature_selector="locus_tag=alpha", placement=side),
        dict(feature_selector="locus_tag=beta", placement="main"),
    ])
    plan = plan_request(_request(mode=mode, feature_placement_table=rows))
    result = _assignments(_planned_layers(plan, mode=mode, direction=direction, resolve=resolve))
    assert (result[1].side, result[1].level) == (side, 1)
    assert result[2].level == 0
    assert len(result) == 3  # default repeat underlay is not a placement unit.


@pytest.mark.parametrize("side", ["outward", "inward"])
@pytest.mark.parametrize("resolve", [False, True])
def test_separate_circular_fixed_placement_uses_the_requested_physical_pool(side, resolve):
    source = _record()
    # Opposite biological strands assigned to one lane still share occupancy.
    source.features[2].location = SimpleLocation(20, 40, -1)
    table = DataFrame([
        dict(feature_selector="locus_tag=alpha", placement=side),
        dict(feature_selector="locus_tag=beta", placement=side),
    ])
    plan = plan_request(_request(source, feature_placement_table=table))
    with pytest.raises(ValidationError, match="Fixed feature placement conflict"):
        _planned_layers(plan, separate=True, resolve=resolve)
    plan = plan_request(_request(source, feature_placement_table=table.iloc[:1]))
    result = _assignments(_planned_layers(plan, separate=True, resolve=resolve))
    assert (result[1].side, result[1].level) == (side, 1)
    assert result[1].strand_pool == ("positive" if side == "outward" else "negative")
    assert result[2].level == 0
    assert result[2].strand_pool == "negative"


@pytest.mark.parametrize("mode,direction,separate", [
    ("circular", "inside", False), ("circular", "outside", False),
    ("circular", "inside", True), ("circular", "outside", True),
    ("linear", "above", False), ("linear", "below", False),
    ("linear", "overlay", True), ("linear", "above", True), ("linear", "below", True),
])
@pytest.mark.parametrize("negative", [False, True])
def test_unsupported_resolved_slot_directions_are_errors(mode, direction, separate, negative):
    side = (("outward", "inward") if mode == "circular" else ("above", "below"))[negative]
    plan = plan_request(_request(mode=mode, feature_placement_table=_table(placement=side)))
    with pytest.raises(ValidationError, match="unsupported"):
        _planned_layers(plan, mode=mode, direction=direction, separate=separate)


@pytest.mark.parametrize("resolve", [False, True])
@pytest.mark.parametrize("placement", ["main", "outward", "inward"])
@pytest.mark.parametrize("tolerance", [19, 20, 21])
def test_fixed_conflict_threshold_is_independent_of_resolver(resolve, placement, tolerance):
    table = DataFrame([
        dict(feature_selector="locus_tag=alpha", placement=placement),
        dict(feature_selector="locus_tag=beta", placement=placement),
    ])
    plan = plan_request(_request(feature_placement_table=table))
    if tolerance < 20:
        with pytest.raises(ValidationError, match="Fixed feature placement conflict"):
            _planned_layers(plan, resolve=resolve, tolerance=tolerance)
    else:
        result = _assignments(_planned_layers(plan, resolve=resolve, tolerance=tolerance))
        assert result[1].level == result[2].level


@pytest.mark.parametrize("resolve", [False, True])
def test_fixed_opposite_sides_can_overlap_and_auto_off_keeps_nominal(resolve):
    plan = plan_request(_request(feature_placement_table=DataFrame([
        dict(feature_selector="locus_tag=alpha", placement="outward"),
        dict(feature_selector="locus_tag=beta", placement="inward"),
    ])))
    result = _assignments(_planned_layers(plan, resolve=resolve))
    assert (result[1].side, result[2].side) == ("outward", "inward")
    assert result[3].level == 0


def test_fresh_auto_removal_and_a_b_a_never_reuse_derived_lanes():
    request = _request(feature_placement_table=_table("locus_tag=beta"))
    first = _planned_layers(plan_request(request))
    beta_main = _assignments(first)
    # Corrupt only a derived runtime carrier to prove neither plan nor next run seeds it.
    next(iter(first.foreground_features.values())).feature_track_id = 92
    auto_plan = plan_request(replace(request, options=replace(
        request.options, feature_placement_table=_table("locus_tag=beta", "auto"),
    )))
    assert auto_plan.inputs.placements == ()
    auto = _assignments(_planned_layers(auto_plan))
    assert (auto[1].level, auto[2].level) == (0, 1)
    last = _planned_layers(plan_request(request))
    assert _assignments(last) == beta_main
    assert all(a is not b for a, b in zip(first.foreground_features.values(), last.foreground_features.values()))


@pytest.mark.parametrize("start", [None, 1, 51])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("crop", [False, True])
def test_transformed_binding_occupancy_and_source_nonmutation(start, reverse, crop):
    source = _record()
    before = copy.deepcopy(source)
    record = RecordInput(
        InMemoryRecordSource(source), record_key="one",
        presentation=RecordPresentation(reverse_complement=reverse),
        region=parse_region_spec("10-105") if crop else None,
        display=RecordDisplayOptions(start_coordinate=start),
    )
    request = _request(records=[record], feature_placements=(_exact(source, 2),))
    if crop and start is not None:
        with pytest.raises(ValidationError, match="explicit display start.*crop"):
            plan_request(request)
        return
    plan = plan_request(request)
    result = _assignments(_planned_layers(plan))
    assert (result[1].level, result[2].level) == (1, 0)
    assert len(result) == 3
    assert str(source.seq) == str(before.seq)
    assert source.annotations == before.annotations
    assert source.features == before.features


def test_nested_multipart_and_artificial_display_seam_are_one_envelope():
    source = _record()
    source.features = [
        SeqFeature(CompoundLocation([SimpleLocation(10, 20, 1), SimpleLocation(90, 110, 1)]),
                   type="CDS", qualifiers={"locus_tag": ["parent"]}),
        SeqFeature(SimpleLocation(40, 60, 1), type="CDS", qualifiers={"locus_tag": ["nested"]}),
    ]
    request = _request(records=[RecordInput(
        InMemoryRecordSource(source), record_key="one", display=RecordDisplayOptions(start_coordinate=16),
    )], feature_placements=(_exact(source, 0),))
    layers = _planned_layers(plan_request(request))
    features = {f.source_feature_index: f for f in layers.foreground_features.values()}
    assert len(features) == 2
    assert len(features[0].display_parts) > 2
    assert (features[0].placement.level, features[1].placement.level) == (0, 1)


@pytest.mark.parametrize("status", ["hidden", "underlay", "crop_excluded"])
def test_dormant_bindings_reserve_no_foreground_lane_and_reactivate(status):
    source = _record()
    target_index = 4 if status == "underlay" else 2
    exact = _exact(source, target_index)
    options = {"feature_placements": (exact,)}
    record = RecordInput(InMemoryRecordSource(source), record_key="one")
    if status == "hidden":
        options["feature_placements"] = (_exact(source, 0),)
        target_index = 0
    elif status == "crop_excluded":
        record = replace(record, region=parse_region_spec("41-110"))
    plan = plan_request(_request(records=[record], **options))
    item = plan.inputs.placements[0].overrides[0]
    assert item.status == status
    layers = _planned_layers(plan)
    assert all(f.placement.requested_target is None for f in layers.foreground_features.values())
    if status != "hidden":
        active = plan_request(_request(feature_placements=(exact,), feature_shapes={"repeat_region": "rectangle"}))
        assert _assignments(_planned_layers(active))[target_index].requested_target == exact.target


@pytest.mark.parametrize("mode,direction", [("circular", "split"), ("linear", "overlay")])
def test_actual_typed_build_transports_instance_context_and_custom_slot(monkeypatch, mode, direction):
    from gbdraw.features import placement
    from gbdraw.api.options import CircularTrackOptions, LinearTrackOptions
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

    source = _record()
    records = [RecordInput(InMemoryRecordSource(source), record_key=key) for key in ("left", "right")]
    slot = (
        CircularTrackSlot(id="custom", renderer="features", params={"lane_direction": direction})
        if mode == "circular" else LinearTrackSlot(id="custom", renderer="features", side=direction)
    )
    track_options = (
        CircularTrackOptions(circular_track_slots=(slot,)) if mode == "circular" else LinearTrackOptions(linear_track_slots=(slot,))
    )
    plan = plan_request(_request(
        mode=mode, records=records, feature_placements=(_exact(source, 2, key="right"),),
        tracks=track_options,
        config_overrides={f"labels.{mode}.scope": "none",
                          "canvas.resolve_overlaps": True, "canvas.strandedness": False},
    ))
    calls = []
    actual = placement.plan_feature_placements

    def capture(features, **kwargs):
        result = actual(features, **kwargs)
        calls.append((kwargs["placement_inputs"].record_key, kwargs["slot"], {
            f.source_feature_index: f.placement for f in features.values()
        }))
        return result

    monkeypatch.setattr(placement, "plan_feature_placements", capture)
    from xml.etree import ElementTree as ET
    built = plan.build()
    root = ET.fromstring((built.drawing if hasattr(built, "drawing") else built).tostring())
    paths = {index: [n.get("d") for n in root.iter()
        if n.get("data-gbdraw-feature-part") == "block"
        and n.get("data-gbdraw-source-feature-index") == index] for index in ("1", "2")}
    assert len(paths["1"]) == len(paths["2"]) == 2
    assert paths["1"][0] == paths["2"][1]
    assert paths["1"][1] == paths["2"][0]
    assert paths["1"][0] != paths["1"][1]
    assert [key for key, _, _ in calls] == ["left", "right"]
    assert [slot.direction for _, slot, _ in calls] == [direction, direction]
    assert [(values[1].level, values[2].level) for _, _, values in calls] == [(0, 1), (1, 0)]


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("form", ["raw", "typed", "overrides"])
def test_nondefault_tolerance_preserves_current_wire_and_copy(mode, form):
    from gbdraw.api.config import load_default_config
    from gbdraw.config.models import GbdrawConfig
    options = {}
    raw = load_default_config()
    raw["canvas"]["feature_overlap_tolerance_bp"] = 2
    if form == "raw":
        options["config"] = raw
    elif form == "typed":
        options["config"] = GbdrawConfig.from_dict(raw)
    else:
        options["config_overrides"] = {"canvas.feature_overlap_tolerance_bp": 2}
    request = _request(mode=mode, **options)
    cloned = replace(request)
    assert cloned.options == request.options
    options_payload = encode_canonical_request(cloned).payload["diagramOptions"]
    assert (options_payload["configOverrides"]["canvas.feature_overlap_tolerance_bp"]
            if form == "overrides" else options_payload["config"]["canvas"]["feature_overlap_tolerance_bp"]) == 2


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_default_tolerance_has_canonical_canvas_owner(mode):
    from gbdraw.api.config import load_default_config
    from gbdraw.config.models import GbdrawConfig
    encoded = encode_canonical_request(_request(mode=mode, config=GbdrawConfig.from_dict(load_default_config())))
    assert encoded.payload["diagramOptions"]["config"]["canvas"]["feature_overlap_tolerance_bp"] == 0


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_manual_origin_wrap_envelope_and_opposite_strand_pool(mode):
    source = _record()
    source.features = [
        SeqFeature(CompoundLocation([SimpleLocation(110, 120, -1), SimpleLocation(0, 10, -1)]),
                   type="CDS", qualifiers={"locus_tag": ["origin"]}),
        SeqFeature(SimpleLocation(2, 8, -1), type="CDS", qualifiers={"locus_tag": ["left"]}),
        SeqFeature(SimpleLocation(2, 8, 1), type="CDS", qualifiers={"locus_tag": ["other"]}),
    ]
    plan = plan_request(_request(source, mode=mode, feature_placements=(_exact(source, 0),)))
    values = _assignments(_planned_layers(
        plan, mode=mode, direction="split" if mode == "circular" else "overlay",
        separate=True,
    ))
    assert [(values[i].strand_pool, values[i].level) for i in range(3)] == [
        ("negative", 0), ("negative", 1), ("positive", 0),
    ]


def test_fixed_first_skips_occupied_directional_lane_and_preserves_auto_priority():
    source = _record()
    source.features = [
        SeqFeature(SimpleLocation(10, 110, 1), type="CDS", qualifiers={"locus_tag": ["long"]}),
        SeqFeature(SimpleLocation(20, 100, 1), type="CDS", qualifiers={"locus_tag": ["auto"]}),
        SeqFeature(SimpleLocation(30, 90, 1), type="CDS", qualifiers={"locus_tag": ["fixed"]}),
    ]
    plan = plan_request(_request(source, feature_placements=(
        _exact(source, 2, target=FeaturePlacementTarget("lane", "outward", 1)),
    )))
    values = _assignments(_planned_layers(plan))
    assert [values[i].level for i in range(3)] == [0, 2, 1]


def test_circular_batch_uses_each_item_placement_context(monkeypatch):
    from gbdraw.features import placement
    source = _record()
    request = CircularBatchRequest(
        records=tuple(RecordInput(InMemoryRecordSource(source), record_key=key)
                      for key in ("left", "right")),
        outputs=tuple(RenderOutputRequest(output_prefix=key) for key in ("left", "right")),
        options=CircularDiagramOptions(
            feature_placements=(_exact(source, 2, key="right"),),
            config_overrides={"canvas.resolve_overlaps": True, "canvas.strandedness": False,
                              "labels.circular.scope": "none", "canvas.circular.track_type": "middle"},
        ),
    )
    calls = []
    actual = placement.plan_feature_placements

    def capture(features, **kwargs):
        result = actual(features, **kwargs)
        calls.append((kwargs["placement_inputs"].record_key,
                      [f.placement.level for f in features.values()]))
        return result

    monkeypatch.setattr(placement, "plan_feature_placements", capture)
    batch = plan_request(request)
    from xml.etree import ElementTree as ET
    paths = []
    for item in batch.item_plans():
        root = ET.fromstring(item.build().tostring())
        paths.append({node.get("data-gbdraw-source-feature-index"): node.get("d")
            for node in root.iter() if node.get("data-gbdraw-feature-part") == "block"
            and node.get("data-gbdraw-source-feature-index") in {"1", "2"}})
    assert paths[0]["1"] == paths[1]["2"]
    assert paths[0]["2"] == paths[1]["1"]
    assert paths[0]["1"] != paths[0]["2"]
    assert calls == [("left", [0, 1, 0]), ("right", [1, 0, 0])]


@pytest.mark.parametrize("mode,side,direction", [
    ("circular", "outward", "split"), ("circular", "inward", "split"),
    ("linear", "above", "overlay"), ("linear", "below", "overlay"),
])
def test_negative_multipart_requested_side_is_independent_of_native_strand(mode, side, direction):
    plan = plan_request(_request(mode=mode, feature_placement_table=_table("locus_tag=multipart", side)))
    layers = _planned_layers(plan, mode=mode, direction=direction)
    feature = next(f for f in layers.foreground_features.values() if f.source_feature_index == 3)
    assert (feature.placement.side, feature.placement.level) == (side, 1)
    assert len(feature.location) == 3  # two blocks and their connector, one assignment.


@pytest.mark.parametrize("mode,direction", [("circular", "inside"), ("linear", "above")])
def test_actual_build_rejects_custom_one_sided_slot_even_under_middle_preset(mode, direction):
    from gbdraw.api.options import CircularTrackOptions, LinearTrackOptions
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot
    tracks = (
        CircularTrackOptions(circular_track_slots=(
            CircularTrackSlot(id="custom", renderer="features", params={"lane_direction": direction}),
        )) if mode == "circular" else LinearTrackOptions(linear_track_slots=(
            LinearTrackSlot(id="custom", renderer="features", side=direction),
        ))
    )
    options = {
        "tracks": tracks,
        "feature_placement_table": _table(placement="outward" if mode == "circular" else "below"),
        "config_overrides": {"canvas.strandedness": False, f"labels.{mode}.scope": "none",
                             "canvas.circular.track_type" if mode == "circular" else "canvas.linear.track_layout": "middle"},
    }
    with pytest.raises(ValidationError, match="unsupported"):
        plan_request(_request(mode=mode, **options)).build()


@pytest.mark.parametrize("form", ["config", "configOverrides"])
def test_historical_reader_does_not_admit_nondefault_tolerance(form, tmp_path):
    from gbdraw.session_request_codec import CanonicalRequestDecodingError, decode_canonical_request
    encoded = encode_canonical_request(_request())
    for resource in encoded.resources:
        (tmp_path / resource.resource_id).write_bytes(resource.content)
    payload = copy.deepcopy(encoded.payload)
    payload["schema"] = 6
    payload["records"][0].pop("display")
    payload["diagramOptions"].pop("featurePlacements")
    payload["diagramOptions"][form] = (
        {"canvas": {"feature_overlap_tolerance_bp": 2}} if form == "config"
        else {"canvas.feature_overlap_tolerance_bp": 2}
    )
    with pytest.raises(CanonicalRequestDecodingError, match="tolerance.*schema 6"):
        decode_canonical_request(
            payload,
            resource_paths={r.resource_id: tmp_path / r.resource_id for r in encoded.resources},
            output_directory=tmp_path,
        )


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_actual_build_requires_an_existing_feature_slot(mode):
    from gbdraw.api.options import CircularTrackOptions, LinearTrackOptions
    from gbdraw.tracks import LinearTrackSlot
    tracks = (
        CircularTrackOptions(circular_track_slots=()) if mode == "circular"
        else LinearTrackOptions(linear_track_slots=(
            LinearTrackSlot(id="gc", renderer="dinucleotide_content", side="below"),
        ))
    )
    plan = plan_request(_request(mode=mode, tracks=tracks, feature_placement_table=_table()))
    with pytest.raises(ValidationError, match="requires resolved slot geometry"):
        plan.build()


@pytest.mark.parametrize("mode,side,style,start,reverse", [
    ("circular", "outward", "horizontal", None, False),
    ("circular", "inward", "horizontal", None, False),
    ("linear", "above", "auto", None, False),
    ("linear", "below", "auto", None, False),
    ("circular", "inward", "radial", 1301, False),
    ("circular", "outward", "radial", 1301, True),
    ("linear", "above", "above_feature", 1301, False),
    ("linear", "below", "above_feature", 1301, True),
])
@pytest.mark.parametrize("resolve", [False, True])
def test_final_svg_isolated_secondary_lane_keeps_empty_main(monkeypatch, mode, side, style, start, reverse, resolve):
    """Measure actual typed assembly and read glyph coordinates from its final SVG."""
    import math
    import re
    from xml.etree import ElementTree as ET
    from gbdraw.api.options import CircularTrackOptions, LinearTrackOptions
    from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot, ScalarSpec
    import gbdraw.diagrams.circular.assemble as circular
    import gbdraw.diagrams.linear.assemble as linear

    record = _record()
    record.seq = Seq("ATGC" * 3000)
    record.features = [SeqFeature(CompoundLocation([
        SimpleLocation(1000, 1800, 1), SimpleLocation(2600, 3200, 1),
    ]), type="CDS", qualifiers={"locus_tag": ["multipart enzyme"]})]
    captures = []
    owner = circular if mode == "circular" else linear
    symbol = "prepare_label_list" if mode == "circular" else "_precalculate_label_dimensions"
    real = getattr(owner, symbol)
    def capture(*args, **kwargs):
        result = real(*args, **kwargs)
        captures.append((kwargs, result))
        return result
    monkeypatch.setattr(owner, symbol, capture)
    if mode == "circular":
        slot = CircularTrackSlot(id="custom", renderer="features", side="overlay",
            radius=ScalarSpec(300, "px"), width=ScalarSpec(20, "px"),
            params={"lane_direction": "split"})
        tracks = CircularTrackOptions(circular_track_slots=(slot,))
    else:
        tracks = LinearTrackOptions(linear_track_slots=(LinearTrackSlot(
            id="custom", renderer="features", side="overlay"),))
    before = copy.deepcopy(record)
    request = _request(record, mode=mode, tracks=tracks,
        records=(RecordInput(InMemoryRecordSource(record), record_key="one",
            display=RecordDisplayOptions(start_coordinate=start),
            presentation=RecordPresentation(reverse_complement=reverse)),),
        feature_placements=(_exact(record, 0, target=FeaturePlacementTarget("lane", side, 1)),),
        config_overrides={"canvas.strandedness": False, "canvas.resolve_overlaps": resolve,
            "canvas.show_gc": False, "canvas.show_skew": False,
            f"labels.{mode}.scope": "both" if mode == "circular" else "all",
            f"labels.{mode}.placement": style,
            "labels.rendering": "auto" if style == "above_feature" else "external_only"})
    drawing = plan_request(request).build()
    root = ET.fromstring((drawing.drawing if hasattr(drawing, "drawing") else drawing).tostring())
    paths = [p for p in root.iter() if p.get("data-gbdraw-feature-part") in {"block", "connector"}
             and not p.get("id", "").endswith("__outline")]
    assert len(paths) == (3 if start is None else 4)
    assert record.seq == before.seq and record.features == before.features
    assert record.annotations == before.annotations
    number = r"[-+0-9.eE]+"
    kwargs, result = captures[-1]
    if mode == "circular":
        layout = kwargs["feature_layout"]
        lane = layout.lanes_by_track_id[1 if side == "outward" else -1]
        # Width 20 and the existing axis-derived gap 3.9 give a 23.9 lane step.
        # An occupied inward lane must not displace the empty Main band.
        main_center = 300
        expected_center = 323.9 if side == "outward" else 276.1
        assert lane.center_px == pytest.approx(expected_center)
        assert layout.primary_band_px.center_px == pytest.approx(main_center)
        assert abs(lane.center_px - main_center) == pytest.approx(23.9)
        for path in paths:
            x, y = map(float, re.search(rf"M\s*({number})[,\s]+({number})", path.get("d")).groups())
            assert lane.inner_px - 1e-6 <= math.hypot(x, y) <= lane.outer_px + 1e-6
        label = result[0]
        assert label["is_inner"] == (side == "inward")
        assert math.hypot(label["feature_middle_x"], label["feature_middle_y"]) == pytest.approx(lane.center_px)
        anchor = (label["feature_anchor_x"], label["feature_anchor_y"])
        assert math.hypot(*anchor) == pytest.approx(lane.inner_px if side == "inward" else lane.outer_px)
        assert any(
            (float(line.get("x2")), float(line.get("y2"))) == pytest.approx(anchor)
            for line in root.iter() if line.tag.endswith("}line")
        )
        assert layout.all_band_px.inner_px <= min(main_center - 10, lane.inner_px)
        assert layout.all_band_px.outer_px >= max(main_center + 10, lane.outer_px)
    else:
        geometry = kwargs["feature_lane_geometries"][0]
        lane = geometry.lanes[0]
        assert (lane.bottom_y < 0) if side == "above" else (lane.top_y > 0)
        for path in paths:
            values = [float(y) for _, y in re.findall(rf"[ML]\s*({number})[,\s]+({number})", path.get("d"))]
            assert all(lane.top_y - 1e-6 <= y <= lane.bottom_y + 1e-6 for y in values)
        label = result[1][0][0]
        assert label["feature_middle_y"] == pytest.approx(lane.middle_y)
        assert (label["middle_y"] < lane.top_y) if side == "above" else (label["middle_y"] > lane.bottom_y)
        if not label["is_embedded"]:
            assert any(
                (float(line.get("x1")), float(line.get("y1"))) == pytest.approx((label["middle"], lane.middle_y))
                for line in root.iter() if line.tag.endswith("}line")
            )
        assert geometry.occupied_band.top_y <= lane.top_y
        assert geometry.occupied_band.bottom_y >= lane.bottom_y


def test_radial_preflight_replans_final_slot_without_stale_assignment(monkeypatch):
    from pathlib import Path
    from gbdraw.features import placement
    import gbdraw.diagrams.circular.assemble as assembly
    import gbdraw.labels.circular as labels

    record = SeqIO.read(Path(__file__).parent / "test_inputs/MjeNMV.gbk", "genbank")
    index = next(i for i, feature in enumerate(record.features) if feature.type == "CDS")
    request = _request(record, feature_placements=(_exact(record, index),),
        config_overrides={"labels.circular.scope": "both", "labels.circular.placement": "radial",
            "labels.rendering": "external_only", "canvas.strandedness": False})
    plans, layouts, candidates = [], [], []
    real_plan, real_labels, real_candidates = placement.plan_feature_placements, assembly.prepare_label_list, labels.build_circular_label_candidates
    def capture_plan(features, **kwargs):
        result = real_plan(features, **kwargs)
        plans.append((kwargs["slot"], result, features))
        return result
    def capture_labels(*args, **kwargs):
        result = real_labels(*args, **kwargs)
        layouts.append((kwargs, result))
        return result
    def capture_candidates(*args, **kwargs):
        result = real_candidates(*args, **kwargs)
        candidates.append(result)
        return result
    monkeypatch.setattr(placement, "plan_feature_placements", capture_plan)
    monkeypatch.setattr(assembly, "prepare_label_list", capture_labels)
    monkeypatch.setattr(labels, "build_circular_label_candidates", capture_candidates)
    first = plan_request(request).build().tostring()
    assert [slot.direction for slot, _, _ in plans] == ["inside", "outside"]
    assert len(candidates) == 1  # Text/source candidates are independent of final lane geometry.
    final = plans[-1][1]
    assert all(f.placement is final[key] for key, f in plans[-1][2].items())
    assert all(a.side == "outward" for a in final.values())
    assert all(a is not plans[0][1][key] for key, a in final.items())
    layout = layouts[-1][0]["feature_layout"]
    assert all(label["feature_center_radius_px"] == pytest.approx(
        layout.lane_for_track_id(label["track_id"]).center_px
    ) for label in layouts[-1][1])
    from xml.etree import ElementTree as ET
    root = ET.fromstring(first)
    endpoints = [(float(line.get("x2")), float(line.get("y2")))
                 for line in root.iter() if line.tag.endswith("}line")]
    assert all(any(point == pytest.approx((label["feature_anchor_x"], label["feature_anchor_y"]))
                   for point in endpoints) for label in layouts[-1][1] if not label["is_embedded"])
    feature_slot = next(slot for slot in layouts[-1][0]["radial_layout"].slots if slot.renderer == "features")
    assert feature_slot.reserved_band_px.inner_px <= layout.all_band_px.inner_px
    assert feature_slot.reserved_band_px.outer_px >= layout.all_band_px.outer_px
    assert max(label.get("required_radius_growth_px", 0) for label in layouts[-1][1]) <= 1e-3
    plans.clear()
    second = plan_request(request).build().tostring()
    assert second == first
    assert [slot.direction for slot, _, _ in plans] == ["inside", "outside"]
    assert plans[-1][1] == final


@pytest.mark.parametrize("mode,side", [("circular", "inward"), ("linear", "below")])
def test_final_svg_placement_a_b_a_reuses_source_cache(tmp_path, monkeypatch, mode, side):
    import gbdraw.api.request_render as owner
    from gbdraw.features import placement
    record = _record()
    record.seq = Seq("ATGC" * 3000)
    path = tmp_path / "source.gb"
    SeqIO.write(record, path, "genbank")
    cache = PreparedBiologicalInputCache()
    identity = PreparedResourceIdentity("source", "placement-resource", path.stat().st_size)
    loads, assignments, svgs, plans = [], [], [], []
    real_load, real_plan = owner.load_gbks, placement.plan_feature_placements
    def load(*args, **kwargs):
        loads.append(1)
        return real_load(*args, **kwargs)
    def capture(*args, **kwargs):
        result = real_plan(*args, **kwargs)
        assignments.append(result)
        return result
    monkeypatch.setattr(owner, "load_gbks", load)
    monkeypatch.setattr(placement, "plan_feature_placements", capture)
    base = _request(mode=mode, records=(RecordInput(GenBankInputSource(path), record_key="one",
        display=RecordDisplayOptions(start_coordinate=51)),), config_overrides={
            "canvas.strandedness": False, "canvas.resolve_overlaps": True,
            ("canvas.circular.track_type" if mode == "circular" else "canvas.linear.track_layout"): "middle",
            f"labels.{mode}.scope": "both" if mode == "circular" else "all"})
    for overrides in [(), (_exact(record, 3, target=FeaturePlacementTarget("lane", side, 1)),), ()]:
        request = replace(base, options=replace(base.options, feature_placements=overrides))
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            plan = plan_request(request)
            plans.append(plan)
            drawing = plan.build()
            svgs.append((drawing.drawing if hasattr(drawing, "drawing") else drawing).tostring())
    assert len(loads) == 1
    assert plans[0].records[0] is plans[1].records[0] is plans[2].records[0]
    assert svgs[0] != svgs[1] and svgs[0] == svgs[2]
    assert assignments[0] == assignments[2] and assignments[0] != assignments[1]


@pytest.mark.parametrize("start", [None, 41])
@pytest.mark.parametrize("feature_bound,style,reverse", [
    (False, "ribbon", False), (True, "ribbon", False),
    (True, "curve", False), (True, "curve", True),
])
def test_final_placement_ribbons_follow_paint_with_four_unit_clearance(start, feature_bound, style, reverse):
    import re
    from xml.etree import ElementTree as ET
    from tests.test_linear_multi_record_comparisons import _comparison
    source = _record()
    before = copy.deepcopy(source)
    comparison = _comparison(1, 0) if reverse else _comparison(0, 1)
    if feature_bound:
        from gbdraw.features.ids import compute_feature_hash
        rows = []
        for query_index, subject_index in (((1, 3), (3, 2)) if reverse else ((3, 1), (2, 3))):
            row = comparison.matches.iloc[0].to_dict()
            for role, index, prefix in (("query", query_index, "q"), ("subject", subject_index, "s")):
                feature = source.features[index]
                row.update({f"{role}_feature_index": index,
                            f"{role}_view_feature_svg_id": compute_feature_hash(feature, record_id=source.id),
                            f"{prefix}start": int(feature.location.start) + 1,
                            f"{prefix}end": int(feature.location.end)})
            rows.append(row)
        comparison = replace(comparison, matches=DataFrame(rows))
    matches = comparison.matches.copy(deep=True)
    records = tuple(RecordInput(InMemoryRecordSource(source), record_key=key,
        display=RecordDisplayOptions(start_coordinate=start)) for key in ("left", "right"))
    base = _request(mode="linear", records=records, linear_comparisons=[comparison], pairwise_match_style=style,
        config_overrides={"canvas.strandedness": False, "canvas.resolve_overlaps": False,
            "canvas.linear.track_layout": "middle", "canvas.show_gc": False, "canvas.show_skew": False,
            "labels.linear.scope": "all", "labels.linear.placement": "above_feature",
            "labels.linear.rotation": 45,
            "objects.features.block_stroke_width.short": 2,
            "objects.features.line_stroke_width.short": 2})
    moved = replace(base, options=replace(base.options, feature_placements=(
        _exact(source, 3, "left", FeaturePlacementTarget("lane", "below", 1)),
        _exact(source, 1, "right", FeaturePlacementTarget("lane", "above", 1)),
    )))
    roots = [ET.fromstring(plan_request(request).build().drawing.tostring()) for request in (base, moved)]
    root = roots[-1]
    parents = {child: parent for parent in root.iter() for child in parent}
    number = r"[-+0-9.eE]+"
    def absolute_y(element, y):
        while element is not None:
            transform = element.get("transform", "")
            assert not any(op in transform for op in ("scale", "matrix", "rotate"))
            y += sum(float(dy) for _, dy in re.findall(rf"translate\(\s*({number})[,\s]+({number})\)", transform))
            element = parents.get(element)
        return y
    bands, feature_bands = {}, {}
    for feature in root.iter():
        if feature.get("data-gbdraw-feature-id") is None or feature.get("d") is None:
            continue
        index = int(feature.get("data-gbdraw-record-index"))
        y_values = [absolute_y(feature, float(y)) for _, y in re.findall(rf"[ML]\s*({number})[,\s]+({number})", feature.get("d"))]
        half_stroke = float(feature.get("stroke-width", "0")) / 2
        top, bottom = min(y_values) - half_stroke, max(y_values) + half_stroke
        prior = bands.get(index, (top, bottom))
        bands[index] = (min(top, prior[0]), max(bottom, prior[1]))
        feature_id = feature.get("data-gbdraw-feature-id")
        prior = feature_bands.get(feature_id, (top, bottom))
        feature_bands[feature_id] = (min(top, prior[0]), max(bottom, prior[1]))
        element = feature
        while element is not None:
            element_id = element.get("id")
            if element_id:
                prior = feature_bands.get(element_id, (top, bottom))
                feature_bands[element_id] = (min(top, prior[0]), max(bottom, prior[1]))
            element = parents.get(element)
    paths = [p for p in root.iter() if p.get("data-pairwise-match-style")]
    assert len(paths) == (2 if feature_bound else 1 if start is None else 3)
    for path in paths:
        endpoints = re.findall(rf"[ML]\s*({number})[,\s]+({number})", path.get("d"))
        query_y, subject_y = [absolute_y(path, float(endpoints[i][1])) for i in (0, 1 if style == "curve" else 2)]
        query_band = feature_bands[path.get("data-query-feature-svg-id")] if feature_bound else bands[0]
        subject_band = feature_bands[path.get("data-subject-feature-svg-id")] if feature_bound else bands[1]
        assert (query_band[0] - query_y if reverse else query_y - query_band[1]) == pytest.approx(4)
        assert (subject_y - subject_band[1] if reverse else subject_band[0] - subject_y) == pytest.approx(4)
    def identity(tree):
        return {p.get("data-gbdraw-match-id") for p in tree.iter() if p.get("data-pairwise-match-style")}
    assert identity(roots[0]) == identity(root)
    assert [p.get("d") for p in roots[0].iter() if p.get("data-pairwise-match-style")] != [p.get("d") for p in paths]
    assert comparison.matches.equals(matches)
    assert source.seq == before.seq and source.features == before.features


def test_circular_existing_cli_keeps_separate_strand_resolver(tmp_path, monkeypatch):
    import re
    from xml.etree import ElementTree as ET
    from gbdraw.circular import circular_main
    import math
    from gbdraw.features import placement
    record = _record()
    path = tmp_path / "source.gb"
    SeqIO.write(record, path, "genbank")
    captured = []
    actual = placement.plan_feature_placements
    def capture(features, **kwargs):
        result = actual(features, **kwargs)
        captured.append((kwargs, {f.source_feature_index: f.placement for f in features.values()}))
        return result
    monkeypatch.setattr(placement, "plan_feature_placements", capture)
    args = ["--gbk", str(path), "--separate_strands", "--resolve_overlaps", "--legend", "none",
            "--format", "svg", "--output", str(tmp_path / "separate")]
    circular_main(args)
    assert captured[-1][0]["resolve_overlaps"] is True
    assert captured[-1][1][1].level == 0 and captured[-1][1][2].level == 1
    root = ET.parse(tmp_path / "separate.svg").getroot()
    radii = {}
    for node in root.iter():
        if node.get("data-gbdraw-feature-part") == "block" and node.get("data-gbdraw-source-feature-index") in {"1", "2"}:
            point = re.search(r"M\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)", node.get("d"))
            radii[int(node.get("data-gbdraw-source-feature-index"))] = math.hypot(*map(float, point.groups()))
    assert radii[1] != radii[2]


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("undefined", [None, 0])
@pytest.mark.parametrize("resolve", [False, True])
def test_final_svg_undefined_uses_negative_geometry_pool(mode, undefined, resolve):
    import math
    import re
    from xml.etree import ElementTree as ET
    record = _record()
    record.seq = Seq("ATGC" * 3000)
    record.features = [SeqFeature(SimpleLocation(1000, 2000, strand), type="CDS",
        qualifiers={"locus_tag": [label]}) for strand, label in ((undefined, "undefined"), (-1, "negative"), (1, "positive"))]
    request = _request(record, mode=mode, feature_placements=(_exact(record, 0),),
        config_overrides={"canvas.strandedness": True, "canvas.resolve_overlaps": resolve,
            ("canvas.circular.track_type" if mode == "circular" else "canvas.linear.track_layout"): "middle",
            f"labels.{mode}.scope": "both" if mode == "circular" else "all"})
    drawing = plan_request(request).build()
    root = ET.fromstring((drawing.drawing if hasattr(drawing, "drawing") else drawing).tostring())
    starts = {}
    for node in root.iter():
        if node.get("data-gbdraw-feature-part") != "block":
            continue
        x, y = map(float, re.search(r"M\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)", node.get("d")).groups())
        starts[node.get("data-gbdraw-feature-id")] = math.hypot(x, y) if mode == "circular" else y
    from gbdraw.features.ids import compute_feature_hash
    # The SVG carrier is not canonical source identity (notably raw strand 0).
    values = [starts[compute_feature_hash(feature, record_id=record.id)] for feature in record.features]
    assert (abs(values[0] - values[1]) > 1) if resolve else (values[0] == pytest.approx(values[1]))
    assert (values[2] > values[0]) if mode == "circular" else (values[2] < values[0])


def test_conflicting_custom_slot_is_rejected_before_radial_reflow():
    from gbdraw.api.options import CircularTrackOptions
    from gbdraw.tracks import CircularTrackSlot
    # A supported directional lane cannot enter the inside-only relocation path:
    # explicit side/direction disagreement is rejected at the existing slot boundary.
    with pytest.raises(ValidationError, match="conflicting side='inside'.*lane_direction='split'"):
        CircularTrackOptions(circular_track_slots=(CircularTrackSlot(
            id="features", renderer="features", side="inside", params={"lane_direction": "split"}),))


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("resolve", [False, True])
def test_final_svg_intron_child_reserves_main_against_whole_parent(mode, resolve):
    import math
    import re
    from xml.etree import ElementTree as ET
    from gbdraw.features.ids import compute_feature_hash

    record = _record()
    record.seq = Seq("ATGC" * 3000)
    record.features = [
        SeqFeature(CompoundLocation([SimpleLocation(1000, 1800, 1), SimpleLocation(2600, 3200, 1)]),
                   type="CDS", qualifiers={"locus_tag": ["parent"]}),
        SeqFeature(SimpleLocation(2000, 2400, 1), type="CDS", qualifiers={"locus_tag": ["intron child"]}),
    ]
    before = copy.deepcopy(record)
    request = _request(record, mode=mode, feature_placements=(_exact(record, 1),),
        config_overrides={"canvas.strandedness": False, "canvas.resolve_overlaps": resolve,
            ("canvas.circular.track_type" if mode == "circular" else "canvas.linear.track_layout"): "middle",
            "canvas.show_gc": False, "canvas.show_skew": False,
            f"labels.{mode}.scope": "both" if mode == "circular" else "all"})
    drawing = plan_request(request).build()
    root = ET.fromstring((drawing.drawing if hasattr(drawing, "drawing") else drawing).tostring())
    rows = []
    for feature in record.features:
        feature_id = compute_feature_hash(feature, record_id=record.id)
        paths = [p for p in root.iter() if p.get("data-gbdraw-feature-id") == feature_id
                 and p.get("data-gbdraw-feature-part") in {"block", "connector"}
                 and not p.get("id", "").endswith("__outline")]
        values = []
        for path in paths:
            pairs = re.findall(r"[ML]\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)", path.get("d"))
            values.extend(math.hypot(float(x), float(y)) if mode == "circular" else float(y) for x, y in pairs)
        rows.append((paths, values))
    assert [len(paths) for paths, _ in rows] == [3, 1]
    parent, child = rows[0][1], rows[1][1]
    if resolve:
        assert min(parent) > max(child) if mode == "circular" else max(parent) < min(child)
    else:
        assert min(parent) == pytest.approx(min(child))
        assert max(parent) == pytest.approx(max(child))
    assert record.seq == before.seq and record.features == before.features
