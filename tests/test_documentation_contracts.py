from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CURRENT_TASK_DOCS = (
    "docs/CLI_Reference.md",
    "docs/PYTHON_API.md",
    "docs/FAQ.md",
    "docs/TUTORIALS/4_Protein_Comparisons.md",
    "docs/TUTORIALS/8_Interactive_SVG_Sessions.md",
)


def test_current_task_docs_delegate_persisted_format_details_to_one_authority() -> None:
    forbidden_details = (
        "__gbdraw_legacy_spacing",
        "schema-2 protein",
        "versions 34–38",
        "request schemas 3–4",
        "ui.layoutPreferences",
        "h_[a-z2-7]{26}",
    )

    for relative_path in CURRENT_TASK_DOCS:
        text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
        assert "SESSION_COMPATIBILITY.md" in text
        for detail in forbidden_details:
            assert detail not in text, f"{relative_path} repeats {detail}"

    authority = (REPO_ROOT / "docs/SESSION_COMPATIBILITY.md").read_text(
        encoding="utf-8"
    )
    for detail in forbidden_details:
        assert detail in authority


def test_python_api_describes_its_four_output_forms() -> None:
    text = (REPO_ROOT / "docs/PYTHON_API.md").read_text(encoding="utf-8")

    assert "provides four output forms" in text
    assert "provides three output methods" not in text
    assert "[Build typed render requests](./TYPED_API.md)" in text


def test_web_maintenance_guide_describes_typed_worker_rendering() -> None:
    text = (REPO_ROOT / "gbdraw/web/CLAUDE.md").read_text(encoding="utf-8")

    assert "canonical schema-5 render request" in text
    assert "workers/diagram-generation-worker.js" in text
    assert "typed request decoding and render_request()" in text
    assert "setup() ~" not in text
    assert "args ~" not in text
