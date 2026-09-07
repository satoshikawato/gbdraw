from __future__ import annotations

import json
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.record_planning import resolve_record_inputs
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.io.regions import parse_region_spec
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.session import build_session_document
from gbdraw.web_support.feature_catalog import build_feature_catalog_item
from gbdraw.web_support.request_render import (
    render_embedded_canonical_web_request,
)


def _colliding_record() -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 20),
        id="collision-record",
        annotations={"molecule_type": "DNA"},
    )
    preceding = SeqFeature(
        SimpleLocation(1, 4, strand=1),
        type="gene",
        qualifiers={"label": ["preceding"]},
    )
    preceding.sub_features = [
        SeqFeature(
            SimpleLocation(2, 3, strand=1),
            type="misc_feature",
            qualifiers={"label": ["nested-preceding"]},
        )
    ]
    record.features = [
        preceding,
        SeqFeature(
            SimpleLocation(10, 19, strand=1),
            type="CDS",
            qualifiers={"label": ["alpha"], "translation": ["MMM"]},
        ),
        SeqFeature(
            SimpleLocation(10, 19, strand=1),
            type="CDS",
            qualifiers={"label": ["beta"], "translation": ["MMM"]},
        ),
    ]
    return record


def _resolve(record_input: RecordInput) -> SeqRecord:
    resolved = resolve_record_inputs(
        (record_input,),
        gff_candidate_features=None,
        gff_keep_all_features=True,
    )
    return resolved.records[0]


def _catalog_identity(
    record: SeqRecord,
) -> tuple[dict[str, str], dict[str, str]]:
    context = build_interactive_svg_context(
        [record],
        selected_features_set=["CDS"],
        linear_rendered_feature_ids=True,
    )
    colors = ("#aa0000", "#0000aa")
    paths = "".join(
        (
            f'<path id="{feature["rendered_feature_svg_id"]}" '
            f'data-gbdraw-feature-id="{feature["rendered_feature_svg_id"]}" '
            f'fill="{color}" d="M 0 0 L 1 1" />'
        )
        for feature, color in zip(context.features, colors, strict=True)
    )
    svg = f'<svg xmlns="http://www.w3.org/2000/svg">{paths}</svg>'
    item = build_feature_catalog_item(
        svg,
        context,
        result_index=0,
        result_name="result.svg",
    )
    identities = {
        str(feature["qualifiers"]["label"][0]): str(
            feature["biologicalFeatureId"]
        )
        for feature in item["biologicalFeatures"]
        if "label" in feature.get("qualifiers", {})
    }
    rendered_identities = {
        str(feature["fillColor"]): str(feature["biologicalFeatureId"])
        for feature in item["features"]
    }

    encoded = json.dumps(item)
    assert "_gbdraw_source_feature_index" not in encoded
    assert "_gbdraw_source_feature_location_parts" not in encoded
    return identities, rendered_identities


def _input(
    record: SeqRecord,
    *,
    region: str | None = None,
    reverse_complement: bool = False,
) -> RecordInput:
    return RecordInput(
        source=InMemoryRecordSource(record),
        record_key="stable-record",
        region=parse_region_spec(region) if region else None,
        presentation=RecordPresentation(
            reverse_complement=reverse_complement,
        ),
    )


def test_crop_preserves_collision_disambiguators_from_source_order() -> None:
    source = _colliding_record()
    full_identities, full_rendered = _catalog_identity(
        _resolve(_input(source))
    )
    cropped_identities, cropped_rendered = _catalog_identity(
        _resolve(_input(source, region="13-30"))
    )

    assert cropped_identities["alpha"] == full_identities["alpha"]
    assert cropped_identities["beta"] == full_identities["beta"]
    assert full_identities["alpha"] != full_identities["beta"]
    assert full_rendered == cropped_rendered == {
        "#aa0000": full_identities["alpha"],
        "#0000aa": full_identities["beta"],
    }
    assert all(
        not hasattr(feature, "_gbdraw_source_feature_index")
        for feature in source.features
    )


def test_reverse_complement_preserves_collision_disambiguators_from_source_order() -> None:
    source = _colliding_record()
    full_identities, full_rendered = _catalog_identity(
        _resolve(_input(source))
    )
    reversed_identities, reversed_rendered = _catalog_identity(
        _resolve(_input(source, reverse_complement=True))
    )

    assert reversed_identities == full_identities
    assert full_identities["alpha"] != full_identities["beta"]
    assert full_rendered == reversed_rendered == {
        "#aa0000": full_identities["alpha"],
        "#0000aa": full_identities["beta"],
    }
    assert all(
        not hasattr(feature, "_gbdraw_source_feature_index")
        for feature in source.features
    )


def test_renderer_emits_unique_catalog_refs_for_exact_feature_collisions(
    tmp_path: Path,
) -> None:
    source = _colliding_record()
    document = build_session_document(
        CircularDiagramRequest(
            records=(_input(source),),
            output=RenderOutputRequest(
                output_prefix="collision",
                formats=("svg",),
            ),
        )
    ).to_dict()

    response = render_embedded_canonical_web_request(
        document["renderRequest"],
        resources=document["resources"],
        workspace=tmp_path / "render",
    )
    item = response["metadata"]["featureCatalog"]["items"][0]
    rendered = item["features"]

    assert len(rendered) == 2
    assert len({feature["svgId"] for feature in rendered}) == 2
    assert len(
        {feature["biologicalFeatureId"] for feature in rendered}
    ) == 2
    for feature in rendered:
        assert (
            f'id="{feature["svgId"]}"'
            in response["results"][0]["content"]
        )


def test_shared_source_catalog_preserves_fixed_preplacement_public_ids():
    # Captured from the Session 6 source bytes before the shared helper existed.
    # Do not derive these expected IDs with the helper under test.
    from gbdraw.features.source import build_source_feature_catalog
    expected = {"preceding": "f86dce137", "nested-preceding": "fe1bb98ed",
                "alpha": "f1ec8c56e~2", "beta": "f1ec8c56e~3"}
    source = _colliding_record()
    catalog = build_source_feature_catalog(source)
    actual = {dict(item.qualifiers)["label"][0]: item.biological_feature_id for item in catalog}
    assert actual == expected
    assert _catalog_identity(_resolve(_input(source)))[0] == expected


def test_gff_collision_ids_use_complete_source_order_across_visibility(tmp_path, monkeypatch):
    from dataclasses import replace

    from gbdraw.api.options import LinearDiagramOptions
    from gbdraw.api.request_render import plan_linear_request
    from gbdraw.api.requests import GffFastaInputSource, LinearDiagramRequest
    from gbdraw.features.placement import FeaturePlacementOverride, FeaturePlacementTarget
    from gbdraw.io import genome

    source = _colliding_record()
    for feature in (*source.features, *source.features[0].sub_features):
        if not hasattr(feature, "sub_features"):
            feature.sub_features = []
    monkeypatch.setattr(genome.GFF, "parse", lambda _: iter((source,)))
    gff, fasta = tmp_path / "source.gff", tmp_path / "source.fa"
    gff.write_text("##gff-version 3\n")
    fasta.write_text(f">{source.id}\n{source.seq}\n")
    record_input = RecordInput(GffFastaInputSource(gff, fasta), record_key="stable-record")
    # The accepted GFF correction replaces filtered-order ~0/~1 with complete
    # source-order ~2/~3. These fixed values come from the preplacement source
    # identity fixture, not from the new shared helper under test.
    expected = {"alpha": "f1ec8c56e~2", "beta": "f1ec8c56e~3"}
    override = FeaturePlacementOverride("stable-record", expected["alpha"], FeaturePlacementTarget("main"))
    for selected in (("CDS",), ("gene", "misc_feature", "CDS")):
        request = LinearDiagramRequest(
            records=(record_input,),
            options=LinearDiagramOptions(selected_features_set=selected),
        )
        plan = plan_linear_request(request)
        ids = _catalog_identity(plan.records[0])[0]
        assert {label: ids[label] for label in expected} == expected
        catalog = plan.provenance[0].source_feature_catalog
        assert [item.biological_feature_id for item in catalog][-2:] == list(expected.values())
        manual = plan_linear_request(replace(request, options=replace(request.options, feature_placements=(override,))))
        bound = manual.inputs.placements[0].foreground
        assert len(bound) == 1
        assert (bound[0].biological_feature_id, bound[0].source_feature_index) == (expected["alpha"], 2)
        from gbdraw.features.factory import create_feature_layers
        from gbdraw.features.placement import FeaturePlacementSlot
        layers = create_feature_layers(
            manual.records[0], {}, ["CDS"], {"default": "#999999", "CDS": "#999999"},
            False, True, {}, compute_label_text=False,
            placement_inputs=manual.inputs.placements[0],
            placement_slot=FeaturePlacementSlot("linear", "overlay", False),
        )
        runtime = {feature.source_feature_index: feature for feature in layers.foreground_features.values()}
        assert (runtime[2].placement.level, runtime[3].placement.level) == (0, 1)
        assert runtime[2].placement.requested_target == override.target

    assert all(not hasattr(feature, "_gbdraw_source_feature_index") for feature in source.features)
