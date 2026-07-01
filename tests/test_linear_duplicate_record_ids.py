import copy
import json
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import gbdraw.linear as linear_cli
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import FeatureDrawingConfigurator
from gbdraw.diagrams.linear.orthogroup_alignment import (
    build_orthogroup_label_eligibility,
    orthogroup_label_sets_for_record,
)
from gbdraw.diagrams.linear.precalc import _precalculate_label_dimensions
from gbdraw.features.ids import compute_feature_hash, make_linear_rendered_feature_id
from gbdraw.io.colors import load_default_colors
from gbdraw.render.interactive_svg import enrich_svg


SVG_NS = "{http://www.w3.org/2000/svg}"


def _record() -> SeqRecord:
    record = SeqRecord(Seq("ATG" * 120), id="duplicate", name="duplicate")
    record.features = [
        SeqFeature(
            SimpleLocation(30, 120, strand=1),
            type="CDS",
            qualifiers={"gene": ["geneA"], "translation": ["M" * 30]},
        ),
        SeqFeature(
            SimpleLocation(180, 270, strand=-1),
            type="CDS",
            qualifiers={"gene": ["geneB"], "translation": ["M" * 30]},
        ),
    ]
    return record


def _reverse_record(record: SeqRecord) -> SeqRecord:
    return record.reverse_complement(id=record.id, name=record.name, description=record.description)


def _depth_table(record_id: str = "duplicate") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [record_id, record_id, record_id],
            "position": [1, 120, 240],
            "depth": [5.0, 20.0, 10.0],
        }
    )


def _linear_feature_config(records: list[SeqRecord]) -> tuple[dict, GbdrawConfig, LinearCanvasConfigurator, FeatureDrawingConfigurator]:
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels="all",
        label_blacklist="",
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=len(records),
        longest_genome=max(len(record.seq) for record in records),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    feature_config = FeatureDrawingConfigurator(
        color_table=None,
        default_colors=load_default_colors("", "default"),
        selected_features_set=cfg.objects.features.features_drawn,
        config_dict=config_dict,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    return config_dict, cfg, canvas_config, feature_config


def _metadata_payload(svg: str) -> dict:
    root = ET.fromstring(svg)
    node = root.find(f".//{SVG_NS}metadata[@id='gbdraw-interactive-feature-metadata']")
    assert node is not None and node.text
    return json.loads(node.text)


def _orthogroup_member(record_index: int, feature_svg_id: str) -> SimpleNamespace:
    return SimpleNamespace(
        orthogroup_id="og_dup",
        protein_id=f"p{record_index}",
        source_protein_id=f"p{record_index}",
        record_index=record_index,
        record_id="duplicate",
        feature_index=0,
        label=f"p{record_index}",
        feature_svg_id=feature_svg_id,
        start=30,
        end=120,
        strand=1,
        representative=record_index == 0,
        role="member",
        confidence="high",
        assignment_reason="test",
        supporting_edges=(),
        best_core_support=0.0,
        second_best_core_support=0.0,
        gene="geneA",
        product=None,
        note=None,
        locus_tag=None,
        gene_id=None,
        old_locus_tag=None,
    )


@pytest.mark.linear
def test_duplicate_record_id_label_precalculation_is_index_aligned_for_reverse_copy() -> None:
    original = _record()
    records = [original, _reverse_record(original)]
    config_dict, cfg, canvas_config, feature_config = _linear_feature_config(records)

    _required_height, labels_by_record, label_heights = _precalculate_label_dimensions(
        records,
        feature_config,
        canvas_config,
        config_dict,
        cfg=cfg,
    )

    assert len(labels_by_record) == 2
    assert len(label_heights) == 2
    assert labels_by_record[0]
    assert labels_by_record[1]
    assert labels_by_record[0] != labels_by_record[1]


@pytest.mark.linear
def test_duplicate_record_id_linear_svg_ids_are_instance_safe() -> None:
    original = _record()
    records = [original, _reverse_record(original)]

    svg = assemble_linear_diagram_from_records(
        records,
        legend="none",
        config_overrides={"show_labels": "all", "show_gc": True, "show_skew": True},
        depth_tables=[_depth_table(), _depth_table()],
        window=60,
        step=60,
        depth_window=60,
        depth_step=60,
    ).tostring()
    root = ET.fromstring(svg)
    ids = [element.attrib["id"] for element in root.iter() if "id" in element.attrib]

    assert len(ids) == len(set(ids))
    assert "duplicate_record_1" in ids
    assert "duplicate_record_2" in ids
    assert "duplicate_definition_record_1" in ids
    assert "duplicate_definition_record_2" in ids
    assert "gc_content_record_1" in ids
    assert "gc_skew_record_2" in ids
    assert "depth_record_1" in ids
    assert "depth_record_2_axis" in ids
    assert svg.count(">duplicate<") == 2


@pytest.mark.linear
def test_duplicate_rendered_feature_payload_entries_do_not_collapse() -> None:
    original = _record()
    records = [original, copy.deepcopy(original)]
    stable_id = compute_feature_hash(original.features[0], record_id=original.id)
    rendered_ids = {
        make_linear_rendered_feature_id(
            record_index=index,
            stable_feature_id=stable_id,
            record_count=len(records),
        )
        for index in range(len(records))
    }

    canvas = assemble_linear_diagram_from_records(
        records,
        legend="none",
        config_overrides={"show_labels": "all"},
    )
    enriched = enrich_svg(
        canvas.tostring(),
        linear_cli._build_interactive_svg_context(records, {"CDS"}),
    )
    payload = _metadata_payload(enriched)
    feature_ids = {feature["svg_id"] for feature in payload["features"]}

    assert rendered_ids.issubset(feature_ids)
    for rendered_id in rendered_ids:
        feature = next(item for item in payload["features"] if item["svg_id"] == rendered_id)
        assert feature["stable_svg_id"] == stable_id


@pytest.mark.linear
def test_orthogroup_label_sets_keep_duplicate_ids_separate_by_record_index() -> None:
    original = _record()
    stable_id = compute_feature_hash(original.features[0], record_id=original.id)
    orthogroups = SimpleNamespace(
        orthogroups={"og_dup": [_orthogroup_member(0, stable_id), _orthogroup_member(1, stable_id)]}
    )

    eligibility = build_orthogroup_label_eligibility(orthogroups=orthogroups)
    lower_member_ids, lower_top_member_ids = orthogroup_label_sets_for_record(eligibility, 1)

    assert stable_id in lower_member_ids
    assert stable_id not in lower_top_member_ids
