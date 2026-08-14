from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import zlib
from pathlib import Path
from xml.etree import ElementTree

import pytest
from PIL import Image

from docs.recipes import run_cli_scenarios as cli_runner
from docs.recipes._scenario_support import (
    PUBLISHED_IMAGE_ROOT,
    extract_executable_block,
    load_chapter,
)
from docs.recipes.run_cli_scenarios import (
    _artifact_comparison_error,
    _assert_png_export,
    _normalized_artifact_payload,
)


pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = "docs/recipes/run_cli_scenarios.py"
SCENARIOS = {
    "H-CLI-02": [
            "record_table.svg",
            "comparison_table.svg",
            "conservation_table.svg",
            "annotation_table.svg",
            "track_table.svg",
    ],
    "H-CLI-09": ["cli_quantitative_tracks.svg"],
    "H-CLI-10": ["cli_annotations_slots.svg"],
    "H-CLI-11": ["cli_feature_presentation.svg"],
    "H-CLI-12": [
        "cli_session.json",
        "cli_session.json.gz",
        "cli_session_roundtrip.svg",
    ],
    "H-CLI-13": [
            "cli_export.svg",
            "cli_export.interactive.svg",
            "cli_export.png",
            "cli_export.pdf",
            "cli_export.eps",
            "cli_export.ps",
    ],
}
EVIDENCE_SOURCE = "docs/internal/SCENARIO_EVIDENCE.md"


def _fixture_manifest() -> dict[str, object]:
    return json.loads(
        (REPO_ROOT / "gbdraw/web/tutorial-data/manifest.json").read_text(
            encoding="utf-8"
        )
    )


def _recipe(scenario_id: str) -> tuple[dict[str, object], str]:
    chapter = load_chapter(
        scenario_id,
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    return chapter, extract_executable_block(chapter, language="bash")


def test_evidence_manifest_entries_are_complete() -> None:
    for scenario_id, expected_outputs in SCENARIOS.items():
        chapter, recipe = _recipe(scenario_id)
        assert chapter["role"] == "evidence"
        assert "destination" not in chapter
        assert chapter["execution"]["source"] == EVIDENCE_SOURCE
        assert chapter["execution"]["expected_outputs"] == expected_outputs
        assert chapter["status"] == {
            "implementation": "verified",
            "review": "approved",
        }
        assert recipe


def test_input_tables_use_whole_sources_and_exact_table_kinds() -> None:
    _, recipe = _recipe("H-CLI-02")
    manifest = _fixture_manifest()
    losatn = manifest["files"]["lambda-de3-losatn"]["expectedSemantics"]
    mt_comparison = manifest["fixtures"]["metazoan-mitochondria-comparison"]

    assert recipe.count("gbdraw ") == 5
    assert "--records_table tables/records.tsv" in recipe
    assert "--comparisons_table tables/comparisons.tsv" in recipe
    assert "--conservation_table tables/conservation.tsv" in recipe
    assert "--annotation_table nicotiana-tabacum-regions.tsv" in recipe
    assert "--circular_track_table tables/tracks.tsv" in recipe
    assert "--gbk NC_001416.gb NC_042057.1.gb" in recipe
    assert "--gbk HmmtDNA.gbk" in recipe
    assert "--conservation_reference subject" in recipe
    assert "--qualifier_priority cds_gene_qualifier_priority.tsv" in recipe
    assert "--labels out" in recipe
    assert "--label_font_size 10" in recipe
    assert manifest["fixtures"]["metazoan-mitochondria-four"][
        "expectedSemantics"
    ]["recordsAreWholeCanonicalSources"]
    assert losatn["queryRecordId"] == "NC_001416.1"
    assert losatn["subjectRecordId"] == "NC_042057.1"
    assert losatn["queryCoordinateRange"] == [1, 48_502]
    assert losatn["subjectCoordinateRange"] == [1, 42_925]
    assert losatn["rawRowCount"] == 6
    assert mt_comparison["expectedSemantics"]["sourceTopologies"] == ["circular"] * 4
    assert mt_comparison["expectedSemantics"]["referenceSide"] == "subject"
    assert mt_comparison["expectedSemantics"]["retainedRows"] == 106


def test_quantitative_and_annotation_recipes_pin_track_semantics() -> None:
    _, depth_recipe = _recipe("H-CLI-09")
    _, annotation_recipe = _recipe("H-CLI-10")
    manifest = _fixture_manifest()

    assert "--depth_window 1" in depth_recipe
    assert "--depth_step 1000" in depth_recipe
    assert "--depth_min 0" in depth_recipe
    assert "--depth_max 80" in depth_recipe
    assert depth_recipe.count("--circular_track_slot") == 6
    assert manifest["fixtures"]["depth-1kb"]["expectedSemantics"][
        "binCount"
    ] == 607

    assert annotation_recipe.count("--circular_track_slot") == 5
    assert "set_id=plastome_regions" in annotation_recipe
    assert manifest["fixtures"]["tobacco-plastome-regions"][
        "expectedSemantics"
    ]["annotationIds"] == ["lsc", "irb", "ssc", "ira"]


def test_presentation_recipe_uses_one_complete_mitochondrion() -> None:
    _, recipe = _recipe("H-CLI-11")

    assert "--palette colorblind" in recipe
    assert "--table tables/presentation_colors.tsv" in recipe
    assert "--label_whitelist tables/presentation_labels.tsv" in recipe
    assert "--label_table tables/presentation_label_overrides.tsv" in recipe
    assert "--feature_visibility_table HmmtDNA_feature_visibility.tsv" in recipe
    assert recipe.count("--feature_shape") == 4
    assert "--resolve_overlaps" in recipe
    assert "--block_stroke_width 1.5" in recipe
    assert "--axis_stroke_width 4" in recipe
    assert "BGC000" not in recipe


def test_session_and_export_recipes_make_overwrite_and_formats_explicit() -> None:
    _, session_recipe = _recipe("H-CLI-12")
    _, export_recipe = _recipe("H-CLI-13")

    assert session_recipe.count("gbdraw circular") == 2
    assert "--session_output cli_session.json" in session_recipe
    assert "--session cli_session.json" in session_recipe
    assert "--session_output cli_session.json.gz" in session_recipe
    assert session_recipe.count("--overwrite") == 1
    assert export_recipe.count("gbdraw circular") == 1
    assert "-f svg,interactive_svg,png,pdf,eps,ps" in export_recipe
    assert SCENARIOS["H-CLI-13"] == [
        "cli_export.svg",
        "cli_export.interactive.svg",
        "cli_export.png",
        "cli_export.pdf",
        "cli_export.eps",
        "cli_export.ps",
    ]


def test_hcli13_fontconfig_uses_default_rules_then_only_bundled_fonts() -> None:
    root = ElementTree.parse(cli_runner._HCLI13_FONTCONFIG_PATH).getroot()
    children = list(root)
    local_names = [element.tag.rsplit("}", 1)[-1] for element in children]

    include = children[local_names.index("include")]
    font_dir = children[local_names.index("dir")]
    cachedir = children[local_names.index("cachedir")]
    assert include.attrib == {"ignore_missing": "no", "prefix": "default"}
    assert include.text == "conf.d"
    assert local_names.index("include") < local_names.index("reset-dirs")
    assert local_names.index("reset-dirs") < local_names.index("dir")
    assert font_dir.attrib == {"prefix": "relative"}
    assert font_dir.text == "../../gbdraw/data"
    assert cachedir.attrib == {"prefix": "xdg"}
    assert cachedir.text == "gbdraw-fontconfig"


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux Fontconfig")
def test_hcli13_fontconfig_resolves_only_bundled_fonts(tmp_path: Path) -> None:
    inherited = os.environ.copy()
    inherited.update(
        {
            "FONTCONFIG_FILE": "/hostile/fonts.conf",
            "FONTCONFIG_PATH": "/hostile/fontconfig",
            "FONTCONFIG_SYSROOT": "/hostile/sysroot",
            "HOME": "/hostile/home",
            "FC_DEBUG": "4096",
            "FC_DBG_MATCH_FILTER": "family",
            "XDG_CACHE_HOME": "/hostile/cache",
            "XDG_CONFIG_HOME": "/hostile/config",
        }
    )
    configured = cli_runner._hcli13_export_environment(
        inherited,
        runtime_root=tmp_path,
    )

    assert inherited["FONTCONFIG_FILE"] == "/hostile/fonts.conf"
    assert configured["FONTCONFIG_FILE"] == str(
        cli_runner._HCLI13_FONTCONFIG_PATH
    )
    assert configured["XDG_CACHE_HOME"] == str(tmp_path / "cache")
    assert configured["XDG_CONFIG_HOME"] == str(tmp_path / "config")
    for name in (
        "FONTCONFIG_PATH",
        "FONTCONFIG_SYSROOT",
        "HOME",
        "FC_DEBUG",
        "FC_DBG_MATCH_FILTER",
    ):
        assert name not in configured

    cli_runner._assert_hcli13_bundled_fonts(configured)
    fc_list = shutil.which("fc-list", path=configured["PATH"])
    assert fc_list is not None
    result = subprocess.run(
        [fc_list, "-f", "%{file}\\n"],
        env=configured,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    resolved = {
        Path(line).resolve() for line in result.stdout.splitlines() if line.strip()
    }
    expected = {path.resolve() for path in cli_runner._HCLI13_FONT_ROOT.glob("*.ttf")}
    assert resolved == expected


def test_hcli13_font_environment_leaves_non_linux_unchanged(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(cli_runner.sys, "platform", "darwin")
    inherited = {
        "HOME": "/original/home",
        "FONTCONFIG_FILE": "/original/fonts.conf",
    }

    configured = cli_runner._hcli13_export_environment(
        inherited,
        runtime_root=tmp_path / "unused",
    )

    assert configured == inherited
    assert configured is not inherited
    assert not (tmp_path / "unused").exists()


def test_hcli13_font_probe_rejects_an_unbundled_resolution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(cli_runner.sys, "platform", "linux")
    monkeypatch.setattr(cli_runner.shutil, "which", lambda *_args, **_kwargs: "fc-match")
    monkeypatch.setattr(
        cli_runner.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                "/host/fonts/LiberationSans-Regular.ttf"
                "\tTrue\tTrue\tFalse\t1\t1\t1\n"
            ),
            stderr="",
        ),
    )

    with pytest.raises(cli_runner.RecipeContractError, match="resolved=/host/fonts"):
        cli_runner._assert_hcli13_bundled_fonts({"PATH": "/host/bin"})


def test_export_png_comparison_tolerates_edge_noise_but_rejects_large_drift(
    tmp_path: Path,
) -> None:
    expected_path = tmp_path / "expected" / "cli_export.png"
    actual_path = tmp_path / "actual" / "cli_export.png"
    expected_path.parent.mkdir()
    actual_path.parent.mkdir()
    expected = Image.new("RGBA", (512, 512), (255, 255, 255, 0))
    expected.paste((40, 100, 180, 255), (40, 40, 472, 472))
    expected.save(expected_path)

    def comparison_error(actual: Image.Image) -> str | None:
        actual.save(actual_path)
        return _artifact_comparison_error("H-CLI-13", expected_path, actual_path)

    assert comparison_error(expected.copy()) is None

    runner_noise = expected.copy()
    for index in range(3_541):
        coordinate = (40 + index % 432, 40 + index // 432)
        runner_noise.putpixel(coordinate, (41, 100, 180, 237))
    runner_noise.putpixel((40, 40), (169, 100, 180, 237))
    assert comparison_error(runner_noise) is None

    recolored = expected.copy()
    recolored.paste((190, 40, 40, 255), (120, 120, 392, 392))
    error = comparison_error(recolored)
    assert error is not None and "large-scale visual content differs" in error

    shifted = Image.new("RGBA", expected.size, (255, 255, 255, 0))
    shifted.paste(expected, (1, 0))
    error = comparison_error(shifted)
    assert error is not None and "large-scale visual content differs" in error

    transparent_blank = Image.new("RGBA", expected.size, (255, 255, 255, 0))
    assert comparison_error(transparent_blank) is not None


def test_assert_png_export_requires_decode_dimensions_mode_and_alpha(
    tmp_path: Path,
) -> None:
    expected_dimensions = (64, 48)
    png_path = tmp_path / "cli_export.png"
    valid = Image.new("RGBA", expected_dimensions, (255, 255, 255, 0))
    valid.paste((20, 80, 160, 255), (8, 8, 56, 40))
    valid.save(png_path)

    _assert_png_export(png_path, expected_dimensions=expected_dimensions)

    truncated_path = tmp_path / "truncated.png"
    truncated_path.write_bytes(png_path.read_bytes()[:-12])
    with pytest.raises(cli_runner.RecipeContractError, match="PNG decode failed"):
        _assert_png_export(
            truncated_path,
            expected_dimensions=expected_dimensions,
        )

    non_png_path = tmp_path / "not-really-png.png"
    valid.save(non_png_path, format="TIFF")
    with pytest.raises(cli_runner.RecipeContractError, match="not a PNG"):
        _assert_png_export(non_png_path, expected_dimensions=expected_dimensions)

    wrong_size_path = tmp_path / "wrong-size.png"
    valid.resize((65, 48)).save(wrong_size_path)
    with pytest.raises(cli_runner.RecipeContractError, match="dimensions"):
        _assert_png_export(wrong_size_path, expected_dimensions=expected_dimensions)

    rgb_path = tmp_path / "rgb.png"
    valid.convert("RGB").save(rgb_path)
    with pytest.raises(cli_runner.RecipeContractError, match="must be RGBA"):
        _assert_png_export(rgb_path, expected_dimensions=expected_dimensions)

    opaque_path = tmp_path / "opaque.png"
    opaque = valid.copy()
    opaque.putalpha(255)
    opaque.save(opaque_path)
    with pytest.raises(cli_runner.RecipeContractError, match="transparent and opaque"):
        _assert_png_export(opaque_path, expected_dimensions=expected_dimensions)


def test_export_renderer_metadata_is_normalized_but_content_is_not(
    tmp_path: Path,
) -> None:
    expected_pdf = tmp_path / "expected.pdf"
    actual_pdf = tmp_path / "actual.pdf"

    def pdf_payload(version: str, date: str, content: str) -> bytes:
        metadata = (
            "1 0 obj\n"
            f"<< /Producer (cairo {version} (https://cairographics.org))\n"
            f"   /CreationDate (D:{date}Z)\n"
            "   /Title (same) >>\n"
            "endobj\n"
        ).encode()
        compressed = zlib.compress(content.encode())
        return (
            metadata
            + b"2 0 obj\n<< /Length 3 0 R >>\nstream\n"
            + compressed
            + b"\nendstream\nendobj\n"
            + f"3 0 obj\n{len(compressed)}\nendobj\n".encode()
        )

    expected_pdf.write_bytes(pdf_payload("1.18.4", "20260804000000", "same"))
    actual_pdf.write_bytes(pdf_payload("1.18.0", "20260810000000", "same"))
    assert _normalized_artifact_payload(
        "H-CLI-13", expected_pdf
    ) == _normalized_artifact_payload("H-CLI-13", actual_pdf)
    assert _artifact_comparison_error(
        "H-CLI-13", expected_pdf, actual_pdf
    ) is None

    actual_pdf.write_bytes(
        actual_pdf.read_bytes().replace(b"cairographics.org", b"example.invalid")
    )
    assert _artifact_comparison_error(
        "H-CLI-13", expected_pdf, actual_pdf
    ) == "normalized payload differs"

    actual_pdf.write_bytes(pdf_payload("1.18.0", "20260810000000", "changed"))
    assert _artifact_comparison_error(
        "H-CLI-13", expected_pdf, actual_pdf
    ) == "normalized payload differs"

    for suffix in (".eps", ".ps"):
        expected_path = tmp_path / f"expected{suffix}"
        actual_path = tmp_path / f"actual{suffix}"
        expected_path.write_bytes(
            b"%%Creator: cairo 1.18.4 (https://cairographics.org)\n"
            b"%%CreationDate: Tue Aug  4 22:21:51 2026\n"
            b"same body\n"
        )
        actual_path.write_bytes(
            b"%%Creator: cairo 1.18.0 (https://cairographics.org)\n"
            b"%%CreationDate: Mon Aug 10 04:43:20 2026\n"
            b"same body\n"
        )
        assert _artifact_comparison_error(
            "H-CLI-13", expected_path, actual_path
        ) is None
        actual_path.write_bytes(
            actual_path.read_bytes().replace(
                b"cairographics.org", b"example.invalid"
            )
        )
        assert _artifact_comparison_error(
            "H-CLI-13", expected_path, actual_path
        ) == "normalized payload differs"
        actual_path.write_bytes(
            actual_path.read_bytes().replace(b"example.invalid", b"cairographics.org")
        )
        actual_path.write_bytes(actual_path.read_bytes().replace(b"same body", b"changed"))
        assert _artifact_comparison_error(
            "H-CLI-13", expected_path, actual_path
        ) == "normalized payload differs"


def test_export_pdf_normalization_respects_declared_stream_length(
    tmp_path: Path,
) -> None:
    expected_pdf = tmp_path / "expected.pdf"
    actual_pdf = tmp_path / "actual.pdf"

    def pdf_payload(version: str, date: str) -> tuple[bytes, bytes]:
        decoded = (
            "16 0\n"
            f"<< /Producer (cairo {version} (https://cairographics.org))\n"
            f"   /CreationDate (D:{date}Z)\n"
            ">>\n"
            "!!!!!!!!!!!!!!!!!!!"
        ).encode()
        compressed = zlib.compress(decoded)
        payload = (
            b"1 0 obj\n<< /Type /ObjStm /Length 2 0 R >>\nstream\n"
            + compressed
            + b"\nendstream\nendobj\n"
            + f"2 0 obj\n{len(compressed)}\nendobj\n".encode()
        )
        return payload, compressed

    expected_payload, expected_stream = pdf_payload("1.18.0", "20260814005730")
    actual_payload, actual_stream = pdf_payload("1.18.4", "20260814005729")
    assert not expected_stream.endswith(b"\r")
    assert actual_stream.endswith(b"\r")
    expected_pdf.write_bytes(expected_payload)
    actual_pdf.write_bytes(actual_payload)

    assert _normalized_artifact_payload(
        "H-CLI-13", expected_pdf
    ) == _normalized_artifact_payload("H-CLI-13", actual_pdf)
    assert _artifact_comparison_error(
        "H-CLI-13", expected_pdf, actual_pdf
    ) is None


@pytest.mark.parametrize("scenario_id", SCENARIOS)
def test_recipe_regenerates_from_a_clean_external_directory(
    scenario_id: str,
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        value
        for value in (str(REPO_ROOT), environment.get("PYTHONPATH"))
        if value
    )
    if scenario_id == "H-CLI-13":
        environment.update(
            {
                "FONTCONFIG_FILE": "/hostile/fonts.conf",
                "FONTCONFIG_PATH": "/hostile/fontconfig",
                "FONTCONFIG_SYSROOT": "/hostile/sysroot",
                "FC_DEBUG": "4096",
                "FC_DBG_MATCH_FILTER": "family",
                "XDG_CACHE_HOME": "/hostile/cache",
                "XDG_CONFIG_HOME": "/hostile/config",
            }
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
        timeout=240,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert f"{scenario_id}: verified" in result.stdout
    assert list(tmp_path.iterdir()) == []
    output_dir = PUBLISHED_IMAGE_ROOT / scenario_id.lower()
    assert all((output_dir / name).is_file() for name in SCENARIOS[scenario_id])
