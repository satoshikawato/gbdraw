"""Collinear option surfaces: one set of defaults and values on every surface (OV-2301..2303)."""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pandas as pd
import pytest

import gbdraw
from gbdraw import LinearComparisonOptions
from gbdraw.analysis.collinearity import (
    COLLINEARITY_ANCHOR_MODES,
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
    LosslessCollinearityParameters,
    convert_collinearity_blocks_to_pair_comparisons,
    normalize_collinearity_anchor_mode,
)
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.exceptions import ValidationError
from gbdraw import linear as linear_cli_module
from gbdraw import session_io
from gbdraw import session_request_codec as codec
from gbdraw.interface import _linear_options

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_collinear_defaults_are_one_set_on_every_surface() -> None:
    """OV-2303: no similarity-group inference and unbounded hit limits by default."""

    assert LinearDiagramOptions().collinear_infer_orthogroups is False
    assert LinearComparisonOptions().collinearity_infer_orthogroups is False
    args = linear_cli_module._get_args(["--gbk", "a.gb", "b.gb", "--losat", "losatp", "--losatp_mode", "collinear"])
    assert args.collinear_infer_orthogroups == "off"
    assert args.losatp_max_target_seqs is None
    assert args.losatp_member_max_hits is None
    assert LosslessCollinearityParameters().merge_orientation == "strand"


def test_beginner_inference_field_reaches_the_typed_options() -> None:
    options = gbdraw.LinearOptions(
        comparisons=LinearComparisonOptions(
            losat="losatp", losatp_mode="collinear", collinearity_infer_orthogroups=True
        )
    )
    assert _linear_options(options, record_count=2).collinear_infer_orthogroups is True


def test_cli_shows_and_wires_every_collinear_setting(capsys: pytest.CaptureFixture[str]) -> None:
    """OV-2303: the three hidden flags are shown; anchor mode and merge orientation are flags."""

    with pytest.raises(SystemExit):
        linear_cli_module._get_args(["--help"])
    help_text = capsys.readouterr().out
    for flag in (
        "--collinear_unit_mode",
        "--collinear_max_conflicts_in_merge_gap",
        "--collinear_max_paralog_links_per_orthogroup",
        "--collinear_anchor_mode {all,rbh}",
        "--collinear_merge_orientation {strand,order}",
    ):
        assert flag in help_text
    assert "color blocks by group" not in help_text

    args = linear_cli_module._get_args([
        "--gbk", "a.gb", "b.gb", "--losat", "losatp", "--losatp_mode", "collinear",
        "--collinear_anchor_mode", "all", "--collinear_merge_orientation", "order",
    ])
    assert (args.collinear_anchor_mode, args.collinear_merge_orientation) == ("all", "order")
    with pytest.raises(SystemExit):
        linear_cli_module._get_args([
            "--gbk", "a.gb", "b.gb", "--losat", "losatp", "--losatp_mode", "collinear",
            "--collinear_anchor_mode", "one_to_one",
        ])


def test_collinear_value_sets_drop_one_to_one_and_read_either_as_strand() -> None:
    """OV-2301, OV-2303: the trimmed value sets."""

    assert COLLINEARITY_ANCHOR_MODES == ("all", "rbh")
    with pytest.raises(ValidationError):
        normalize_collinearity_anchor_mode("one_to_one")
    assert LosslessCollinearityParameters(merge_orientation="either").merge_orientation == "strand"  # type: ignore[arg-type]
    with pytest.raises(ValidationError):
        LosslessCollinearityParameters(merge_orientation="both").validate()  # type: ignore[arg-type]


def _protein_settings(**changes: object) -> dict[str, object]:
    settings = codec._encode_protein_settings(LinearDiagramOptions())
    settings.update(changes)
    return settings


def _decode_settings(settings: dict[str, object]) -> dict[str, object]:
    return codec._decode_pipeline(
        {"kind": "generatedProteinComparison", "mode": "none", "pairs": [], "settings": settings},
        path="comparisons[0]",
        schema=codec.CANONICAL_REQUEST_SCHEMA,
    )


def test_stored_request_reads_retired_values_and_old_inference(caplog: pytest.LogCaptureFixture) -> None:
    """A saved one_to_one reads as rbh with a warning; either as strand; no flag as inference on."""

    params = {
        "kind": "lossless",
        "parameters": {
            "minAnchors": 1, "maxUnitGap": 0, "maxDiagonalDrift": 0, "maxConflicts": 1,
            "mergeOrientation": "either",
        },
    }
    settings = _protein_settings(collinearityAnchorMode="one_to_one", collinearityParams=params)
    settings.pop("collinearInferOrthogroups")
    with caplog.at_level("WARNING"):
        decoded = _decode_settings(settings)
    assert decoded["collinearity_anchor_mode"] == "rbh"
    assert "one_to_one" in caplog.text
    assert decoded["collinearity_params"].merge_orientation == "strand"
    assert decoded["collinear_infer_orthogroups"] is True

    shapes = codec.describe_canonical_request()["types"]
    protein = shapes["CanonicalProteinComparisonSettings"]
    assert protein["defaults"]["collinearInferOrthogroups"] is True
    assert protein["defaults"]["collinearityAnchorMode"] == "rbh"
    merge = shapes["CanonicalLosslessCollinearityParameters"]["fields"]["mergeOrientation"]
    # A released writer stored "either"; the stored value set keeps it, and
    # every reader reads it as "strand".
    assert merge["values"] == ["strand", "order", "either"]
    assert merge["spellings"]["aliases"] == {"either": "strand"}


def _anchor(index: int, *, orthogroup_id: str = "") -> CollinearityAnchor:
    return CollinearityAnchor(
        query_protein_id=f"q{index}", subject_protein_id=f"s{index}",
        query_record_index=0, subject_record_index=1,
        query_order=index, subject_order=index,
        query_start=index * 100 + 1, query_end=index * 100 + 90,
        subject_start=index * 100 + 1, subject_end=index * 100 + 90,
        identity=90.0, evalue=1e-30, bitscore=100.0, alignment_length=30,
        query_feature_svg_id=f"qf{index}", subject_feature_svg_id=f"sf{index}",
        source="rbh", query_unit_id=f"qu{index}", subject_unit_id=f"su{index}",
        query_unit_kind="cds", subject_unit_kind="cds",
        query_locus_id=None, subject_locus_id=None,
        query_display_name=f"q{index}", subject_display_name=f"s{index}",
        query_strand=1, subject_strand=1, orthogroup_id=orthogroup_id,
    )


@pytest.mark.parametrize(("orthogroup_id", "expected_scope"), [("", None), ("og_1", "global_collinear")])
def test_block_rows_name_a_group_scope_only_when_groups_exist(
    orthogroup_id: str, expected_scope: str | None
) -> None:
    """OV-2302: blocks built without similarity groups claim no group scope."""

    block = CollinearityBlock(
        block_id="b1", query_record_index=0, subject_record_index=1, orientation="plus",
        score=200.0, anchors=(_anchor(0, orthogroup_id=orthogroup_id), _anchor(1, orthogroup_id=orthogroup_id)),
        kind="cluster",
    )
    frames = convert_collinearity_blocks_to_pair_comparisons(
        CollinearityResult(blocks=(block,)), record_ids=["A", "B"], search_scope="all"
    )
    frame = frames[(0, 1)]
    if expected_scope is None:
        assert "group_scope" not in frame.columns and "group_kind" not in frame.columns
    else:
        assert set(frame["group_scope"]) == {expected_scope}
        assert set(frame["group_kind"]) == {"orthogroup"}
    assert isinstance(frame, pd.DataFrame)


def _gui_collinear_session(blastp: dict[str, object]) -> dict[str, object]:
    return {"config": {"losat": {"blastp": {"mode": "collinear", **blastp}}}, "version": 30}


@pytest.mark.parametrize(
    ("blastp", "expected"),
    [
        ({}, ["--collinear_infer_orthogroups", "on"]),
        ({"collinearInferOrthogroups": False}, ["--collinear_infer_orthogroups", "off"]),
        (
            {"collinearAnchorMode": "one_to_one", "collinearMergeOrientation": "either"},
            ["--collinear_infer_orthogroups", "on", "--collinear_anchor_mode", "rbh",
             "--collinear_merge_orientation", "strand"],
        ),
    ],
)
def test_legacy_gui_session_replay_keeps_its_collinear_settings(
    blastp: dict[str, object], expected: list[str]
) -> None:
    session = _gui_collinear_session(blastp)
    run_args: list[str] = []
    session_io._append_linear_gui_blastp_args(
        run_args, [], session=session, config=session["config"], adv={},  # type: ignore[arg-type]
    )
    start = run_args.index("--collinear_infer_orthogroups")
    assert run_args[start:] == expected


def test_recorded_collinear_command_without_the_flag_replays_with_inference() -> None:
    args = ["linear", "--gbk", "a.gb", "--losat", "losatp", "--losatp_mode", "collinear"]
    assert session_io._keep_recorded_collinear_inference(args)[-2:] == ["--collinear_infer_orthogroups", "on"]
    explicit = [*args, "--collinear_infer_orthogroups", "off"]
    assert session_io._keep_recorded_collinear_inference(explicit) == explicit
    pairwise = ["linear", "--losat", "losatp", "--losatp_mode", "pairwise"]
    assert session_io._keep_recorded_collinear_inference(pairwise) == pairwise


def _bash_blocks(path: Path) -> list[str]:
    return re.findall(r"```bash\n(.*?)```", path.read_text(encoding="utf-8"), flags=re.S)


def test_every_recorded_collinear_command_states_the_inference_flag() -> None:
    """A published Collinear command reproduces its figure whatever the CLI default."""

    from tools import prepare_interactive_gallery_assets as gallery
    from tools.reproduce_examples_manifest import CliRecipe, build_figure_specs

    commands: dict[str, list[str]] = {}
    for path in [
        *sorted((REPO_ROOT / "docs").rglob("*.md")),
        REPO_ROOT / "docs/internal/SCENARIO_EVIDENCE.md",
    ]:
        if "docs/internal/" in path.as_posix() and path.name != "SCENARIO_EVIDENCE.md":
            continue
        for index, block in enumerate(_bash_blocks(path)):
            tokens = shlex.split(block.replace("\\\n", " "), comments=True)
            if "collinear" in tokens and "--losatp_mode" in tokens:
                commands[f"{path.relative_to(REPO_ROOT)}#{index}"] = tokens
    for example in gallery.EXAMPLES:
        tokens = shlex.split(example.command)
        if "--losatp_mode" in tokens and "collinear" in tokens:
            commands[f"gallery:{example.id}"] = tokens
    for figure_id, spec in build_figure_specs().items():
        if isinstance(spec.recipe, CliRecipe) and "collinear" in spec.recipe.extra_args:
            commands[f"examples:{figure_id}"] = list(spec.recipe.extra_args)

    assert {"gallery:hepatoplasmataceae_collinear", "gallery:vibrio-harveyi-group-collinear"} <= set(commands)
    assert any("compare-proteins-losatp-collinear.md" in name for name in commands)
    assert any("SCENARIO_EVIDENCE.md" in name for name in commands)
    missing = sorted(name for name, tokens in commands.items() if "--collinear_infer_orthogroups" not in tokens)
    assert missing == []


def test_a_released_session_with_the_removed_anchor_mode_loads_as_rbh(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A released schema 9 request may hold one_to_one (OV-2303) and holds "either"."""

    import json

    from gbdraw.session import load_session_document, materialize_session, session_to_request

    path = REPO_ROOT / "tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json"
    session = json.loads(path.read_text(encoding="utf-8"))
    (settings,) = [item["settings"] for item in session["renderRequest"]["comparisons"] if "settings" in item]
    assert settings["collinearityParams"]["parameters"]["mergeOrientation"] == "either"
    settings["collinearityAnchorMode"] = "one_to_one"

    with caplog.at_level("WARNING"), materialize_session(
        load_session_document(session), output_directory=tmp_path
    ) as materialized:
        options = session_to_request(materialized).options
    assert options.collinearity_anchor_mode == "rbh"
    assert options.collinearity_params.merge_orientation == "strand"
    assert "one_to_one" in caplog.text
