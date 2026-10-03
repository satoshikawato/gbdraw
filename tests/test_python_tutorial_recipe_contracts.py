from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
from pandas import DataFrame

import gbdraw.api.diagram as api_diagram_module
from docs.recipes._scenario_support import (
    copy_declared_inputs,
    extract_executable_block,
    load_chapter,
)
from docs.recipes.run_python_scenarios import SCENARIO_IDS
from gbdraw.analysis.protein_colinearity import (
    OrthogroupMember,
    OrthogroupResult,
    ProteinBlastpResult,
)
from gbdraw.api import RequestRenderResult, SimilarityAlignmentPlan
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.layout.similarity_alignment import AlignmentDecisionStatus


pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = "docs/recipes/run_python_scenarios.py"
TUTORIAL_SCENARIO_IDS = tuple(
    scenario_id
    for scenario_id in SCENARIO_IDS
    # Full-data LOSATP recipes remain available through the manual scenario
    # runner (removed from CI for runtime on 2026-08-11). T-PY-05's documented
    # call is covered below with fixed Similarity groups instead of LOSATP.
    if scenario_id.startswith("T-PY-")
    and scenario_id not in {"T-PY-01", "T-PY-05", "T-PY-07"}
)


@pytest.mark.parametrize(
    "scenario_id",
    TUTORIAL_SCENARIO_IDS,
)
def test_python_tutorial_recipe_regenerates_from_a_clean_external_context(
    scenario_id: str,
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not existing_pythonpath
        else os.pathsep.join((str(REPO_ROOT), existing_pythonpath))
    )

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / RUNNER),
            "--scenario",
            scenario_id,
            "--check",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert f"{scenario_id}: verified" in result.stdout
    assert list(tmp_path.iterdir()) == []


def test_t_py_05_program_aligns_on_its_protein_id_without_losatp(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    chapter = load_chapter("T-PY-05", expected_kind="python-recipe", runner_path=RUNNER)
    program = extract_executable_block(chapter, language="python")
    copy_declared_inputs(chapter, recipe_source=program, workdir=tmp_path)
    searches: list[int] = []

    def og_1_with_the_reference_only(records, *, protein_extraction, **_kwargs):
        searches.append(len(records))
        members = [
            OrthogroupMember(
                orthogroup_id="og_1",
                protein_id=protein.protein_id,
                record_index=protein.record_index,
                feature_index=protein.feature_index,
                record_id=protein.record_id,
                label=protein.label,
                start=protein.start,
                end=protein.end,
                strand=protein.strand,
                feature_svg_id=protein.feature_svg_id,
                source_protein_id=protein.source_protein_id,
            )
            for proteins in protein_extraction.proteins_by_record
            for protein in proteins
            if protein.source_protein_id == "CAG38695.1"
        ]
        return ProteinBlastpResult(
            comparisons=[
                DataFrame(columns=COMPARISON_COLUMNS) for _ in range(len(records) - 1)
            ],
            orthogroups=OrthogroupResult(
                orthogroups={"og_1": members},
                member_by_protein_id={member.protein_id: member for member in members},
            ),
        )

    monkeypatch.setattr(
        api_diagram_module,
        "build_rbh_orthogroup_protein_blastp_comparisons",
        og_1_with_the_reference_only,
    )
    monkeypatch.chdir(tmp_path)
    namespace: dict[str, object] = {"__name__": "__gbdraw_documented_recipe__"}
    exec(compile(program, "bgc_losatp_groups.py", "exec"), namespace)

    result = namespace["diagram"]
    assert isinstance(result, RequestRenderResult)
    assert namespace["saved_path"] == Path("python_bgc_losatp_groups.svg")
    assert (tmp_path / "python_bgc_losatp_groups.svg").is_file()
    assert searches == [5]
    plan = result.request.similarity_alignment
    assert isinstance(plan, SimilarityAlignmentPlan)
    assert plan.group_id == "og_1"
    assert plan.reference.record_key == "record-1"
    assert [decision.status for decision in plan.records] == [
        AlignmentDecisionStatus.REFERENCE,
        *([AlignmentDecisionStatus.SKIPPED] * 4),
    ]
    assert result.request.options.losat_search.losatp_mode == "none"
