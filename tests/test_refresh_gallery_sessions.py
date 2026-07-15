from __future__ import annotations

from pathlib import Path

from tools.refresh_gallery_sessions import (
    _prune_derivable_canonical_gallery_fields,
    _preserve_gallery_cli_invocation,
    _with_interactive_svg_format,
)


def test_with_interactive_svg_format_replaces_existing_format() -> None:
    assert _with_interactive_svg_format(["-o", "out", "-f", "svg"]) == [
        "-o",
        "out",
        "-f",
        "interactive_svg",
    ]
    assert _with_interactive_svg_format(["--format=svg"]) == [
        "--format=interactive_svg"
    ]
    assert _with_interactive_svg_format(["--gbk", "input.gb"]) == [
        "--gbk",
        "input.gb",
        "-f",
        "interactive_svg",
    ]


def test_preserve_gallery_cli_invocation_keeps_original_render_args() -> None:
    source_session = {
        "cliInvocation": {
            "schema": 1,
            "mode": "circular",
            "args": [
                "--definition_line_style",
                "name:font_weight=bold",
                "--circular_track_slot",
                "a_skew_2:dinucleotide_skew@w=0.1,nt=AT,legend_label=AT skew",
                "--conservation_blast",
                "comparison.tsv",
                "--gbk",
                "input.gb",
                "-f",
                "svg",
            ],
            "renderFormats": ["svg"],
            "fileBindings": [
                {
                    "argIndex": 5,
                    "slot": "files.c_conservation_blasts[0]",
                    "name": "comparison.tsv",
                },
                {"argIndex": 7, "slot": "files.c_gb", "name": "input.gb"},
            ],
            "generatedBy": "gbdraw",
        },
    }
    refreshed_session = {
        "cliInvocation": {
            "schema": 1,
            "mode": "circular",
            "args": ["--gbk", "input.gb", "-f", "interactive_svg"],
            "renderFormats": ["interactive_svg"],
            "fileBindings": [
                {"argIndex": 1, "slot": "files.c_gb", "name": "input.gb"}
            ],
            "generatedBy": "gbdraw",
        },
    }

    preserved = _preserve_gallery_cli_invocation(
        source_session,
        refreshed_session,
        mode="circular",
    )

    assert preserved is True
    args = refreshed_session["cliInvocation"]["args"]
    assert "--definition_line_style" in args
    assert "--circular_track_slot" in args
    assert any("AT skew" in arg for arg in args)
    assert "--conservation_blast" in args
    assert args[-2:] == ["-f", "interactive_svg"]
    assert refreshed_session["cliInvocation"]["renderFormats"] == ["interactive_svg"]
    assert (
        refreshed_session["cliInvocation"]["fileBindings"]
        == source_session["cliInvocation"]["fileBindings"]
    )


def test_preserve_gallery_cli_invocation_reports_missing_source_cli() -> None:
    refreshed_session = {
        "cliInvocation": {
            "schema": 1,
            "mode": "circular",
            "args": ["--gbk", "input.gb", "-f", "interactive_svg"],
        },
    }

    preserved = _preserve_gallery_cli_invocation(
        {"config": {"form": {"prefix": "old"}}},
        refreshed_session,
        mode="circular",
    )

    assert preserved is False
    assert refreshed_session["cliInvocation"]["args"] == [
        "--gbk",
        "input.gb",
        "-f",
        "interactive_svg",
    ]


def test_prune_derivable_canonical_gallery_fields_keeps_authority_and_preview() -> None:
    session = {
        "version": 31,
        "renderRequest": {"schema": 1},
        "resources": {"record-1-genbank": {}},
        "results": [{"name": "preview", "content": "<svg/>"}],
        "files": {"c_gb": {}},
        "features": {"extractedFeatures": []},
        "losatCache": {"entries": []},
        "losatDerivedCache": {"entries": []},
    }

    _prune_derivable_canonical_gallery_fields(session)

    assert session == {
        "version": 31,
        "renderRequest": {"schema": 1},
        "resources": {"record-1-genbank": {}},
        "results": [{"name": "preview", "content": "<svg/>"}],
    }


def test_prepare_gallery_assets_preserves_existing_source_svgs() -> None:
    source = Path("tools/prepare_interactive_gallery_assets.py").read_text(
        encoding="utf-8"
    )

    assert "def _read_or_create_source_svg(" in source
    assert "if example.source_svg_path.exists():" in source
    assert "example.source_svg_path.write_text(source" in source
    assert "def _sync_session_result_svg(" in source
    assert "write_session_json(example.session_path, session)" in source
    assert "_sync_session_result_svg(example, session, source)" in source
    assert "_write_gallery_svg(example, session, source)" in source
    assert 'entry["tutorial"] = f"./tutorials/{example.id}.json"' in source
    assert 'entry["tutorialStatus"] = "ready"' in source
