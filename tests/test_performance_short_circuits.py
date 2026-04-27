from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.cli_utils.common as cli_common_module
import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.diagrams.linear.assemble as linear_assemble_module
import gbdraw.diagrams.linear.precalc as linear_precalc_module
import gbdraw.features.factory as factory_module
import gbdraw.render.groups.linear.gc_content as linear_gc_content_module
import gbdraw.render.groups.linear.gc_skew as linear_gc_skew_module
import gbdraw.render.groups.linear.seq_record as linear_seq_record_module
import gbdraw.render.export as export_module
from gbdraw.analysis.skew import skew_df as real_skew_df
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)


def _make_record(record_id: str = "rec1") -> SeqRecord:
    record = SeqRecord(Seq("ATGCGT" * 120), id=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(10, 180, strand=1),
            type="CDS",
            qualifiers={
                "gene": ["geneA"],
                "product": ["enzyme alpha"],
            },
        ),
        SeqFeature(
            FeatureLocation(220, 310, strand=-1),
            type="misc_feature",
            qualifiers={"note": ["test feature"]},
        ),
    ]
    return record


@pytest.mark.linear
def test_linear_shared_gc_dataframe_is_computed_once_per_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _make_record()
    call_count = 0

    def counting_skew_df(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return real_skew_df(*args, **kwargs)

    def unexpected_group_skew_df(*_args, **_kwargs):
        raise AssertionError("Linear GC groups should reuse the precomputed GC dataframe.")

    monkeypatch.setattr(linear_assemble_module, "skew_df", counting_skew_df)
    monkeypatch.setattr(linear_gc_content_module, "skew_df", unexpected_group_skew_df)
    monkeypatch.setattr(linear_gc_skew_module, "skew_df", unexpected_group_skew_df)

    assemble_linear_diagram_from_records(
        [record],
        legend="none",
        config_overrides={
            "show_gc": True,
            "show_skew": True,
            "show_labels": "none",
        },
    )

    assert call_count == 1


@pytest.mark.circular
def test_circular_no_labels_skips_label_text_lookup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _make_record()
    call_count = 0

    def counting_get_label_text(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return ""

    monkeypatch.setattr(factory_module, "get_label_text", counting_get_label_text)

    assemble_circular_diagram_from_record(
        record,
        legend="none",
        config_overrides={
            "show_labels": "none",
            "show_gc": False,
            "show_skew": False,
        },
    )

    assert call_count == 0


@pytest.mark.linear
def test_linear_no_labels_skips_label_text_lookup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _make_record()
    call_count = 0

    def counting_get_label_text(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return ""

    monkeypatch.setattr(factory_module, "get_label_text", counting_get_label_text)

    assemble_linear_diagram_from_records(
        [record],
        legend="none",
        config_overrides={
            "show_labels": "none",
            "show_gc": False,
            "show_skew": False,
        },
    )

    assert call_count == 0


@pytest.mark.linear
def test_linear_feature_dict_is_built_once_per_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [_make_record("rec1"), _make_record("rec2")]
    call_count = 0
    real_create_feature_dict = linear_precalc_module.create_feature_dict

    def counting_create_feature_dict(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return real_create_feature_dict(*args, **kwargs)

    def unexpected_create_feature_dict(*_args, **_kwargs):
        raise AssertionError("Linear assembly should reuse precomputed feature dictionaries.")

    monkeypatch.setattr(linear_precalc_module, "create_feature_dict", counting_create_feature_dict)
    monkeypatch.setattr(linear_assemble_module, "create_feature_dict", unexpected_create_feature_dict)
    monkeypatch.setattr(linear_seq_record_module, "create_feature_dict", unexpected_create_feature_dict)

    assemble_linear_diagram_from_records(
        records,
        legend="none",
        config_overrides={
            "show_labels": "all",
            "show_gc": False,
            "show_skew": False,
        },
    )

    assert call_count == len(records)


@pytest.mark.circular
def test_circular_legend_none_skips_legend_preparation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _make_record()

    def unexpected_legend_work(*_args, **_kwargs):
        raise AssertionError("Legend preparation should be skipped when legend='none'.")

    monkeypatch.setattr(circular_assemble_module, "check_feature_presence", unexpected_legend_work)
    monkeypatch.setattr(circular_assemble_module, "precompute_used_color_rules", unexpected_legend_work)
    monkeypatch.setattr(circular_assemble_module, "prepare_legend_table", unexpected_legend_work)

    assemble_circular_diagram_from_record(
        record,
        legend="none",
        config_overrides={
            "show_labels": "none",
            "show_gc": False,
            "show_skew": False,
        },
    )


@pytest.mark.linear
def test_linear_legend_none_skips_legend_preparation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _make_record()

    def unexpected_legend_work(*_args, **_kwargs):
        raise AssertionError("Legend preparation should be skipped when legend='none'.")

    monkeypatch.setattr(linear_assemble_module, "check_feature_presence", unexpected_legend_work)
    monkeypatch.setattr(linear_assemble_module, "precompute_used_color_rules", unexpected_legend_work)
    monkeypatch.setattr(linear_assemble_module, "prepare_legend_table", unexpected_legend_work)

    assemble_linear_diagram_from_records(
        [record],
        legend="none",
        config_overrides={
            "show_labels": "none",
            "show_gc": False,
            "show_skew": False,
        },
    )


def test_save_figure_svg_only_skips_cairosvg_lookup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    drawing = Drawing(filename=str(tmp_path / "figure.svg"))
    drawing.add(drawing.rect(insert=(0, 0), size=("10px", "10px")))

    def unexpected_get_cairosvg():
        raise AssertionError("SVG-only export should not resolve CairoSVG.")

    monkeypatch.setattr(export_module, "get_cairosvg", unexpected_get_cairosvg)

    export_module.save_figure(drawing, ["svg"])

    assert (tmp_path / "figure.svg").exists()


def test_handle_output_formats_svg_only_skips_cairosvg_lookup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def unexpected_has_cairosvg():
        raise AssertionError("SVG-only CLI output should not check CairoSVG availability.")

    monkeypatch.setattr(cli_common_module, "has_cairosvg", unexpected_has_cairosvg)

    assert cli_common_module.handle_output_formats(["svg"]) == ["svg"]
