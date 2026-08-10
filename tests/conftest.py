"""Shared pytest fixtures for gbdraw."""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Callable
from functools import cache
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from gbdraw.session_io import load_session


PROJECT_ROOT = Path(__file__).parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "examples"
TEST_INPUTS_DIR = PROJECT_ROOT / "tests" / "test_inputs"
INPUT_SEARCH_DIRS = (
    TEST_INPUTS_DIR,
    EXAMPLES_DIR,
    Path("/home/kawato/study/2025-09-17_gbdraw_test"),
    Path("/home/kawato/study/2025-05-15_gbdraw"),
)


def pytest_addoption(parser: pytest.Parser) -> None:
    """Register explicit opt-in options for write operations."""

    parser.addoption(
        "--update-reference-outputs",
        action="store_true",
        default=False,
        help="Regenerate tracked SVG reference outputs.",
    )


def pytest_collection_modifyitems(
    config: pytest.Config,
    items: list[pytest.Item],
) -> None:
    """Skip reference generation unless its write operation was requested."""

    if config.getoption("--update-reference-outputs"):
        return

    skip_generation = pytest.mark.skip(
        reason="reference generation requires --update-reference-outputs"
    )
    for item in items:
        if "reference_generation" in item.keywords:
            item.add_marker(skip_generation)


@pytest.fixture(scope="session")
def examples_dir() -> Path:
    return EXAMPLES_DIR


@pytest.fixture(scope="session")
def test_inputs_dir() -> Path:
    return TEST_INPUTS_DIR


@pytest.fixture(scope="session")
def load_cached_gallery_session() -> Callable[[Path], dict[str, object]]:
    @cache
    def load(path: Path) -> dict[str, object]:
        return load_session(path)

    return lambda path: load(path.resolve())


@pytest.fixture(scope="session")
def load_cached_svg_root() -> Callable[[Path], ET.Element]:
    @cache
    def load(path: Path) -> ET.Element:
        return ET.parse(path).getroot()

    return lambda path: load(path.resolve())


@pytest.fixture(scope="session")
def find_test_input() -> Callable[[str], Path | None]:
    """Find an input in the repository or optional local fixture directories."""

    def find(filename: str) -> Path | None:
        return next(
            (
                candidate
                for directory in INPUT_SEARCH_DIRS
                if directory.exists()
                if (candidate := directory / filename).exists()
            ),
            None,
        )

    return find


class GbdrawRunner:
    """Run one gbdraw CLI subcommand and report its SVG output."""

    def run(
        self,
        subcommand: str,
        gbk_files: list[Path],
        output_prefix: str,
        output_dir: Path,
        *,
        blast_files: list[Path] | None = None,
        extra_args: list[str] | None = None,
        timeout: int = 300,
    ) -> tuple[int, str, Path]:
        cmd = [
            sys.executable,
            "-m",
            "gbdraw.cli",
            subcommand,
            "--gbk",
            *(str(path) for path in gbk_files),
            "-o",
            str(output_dir / output_prefix),
            "-f",
            "svg",
        ]
        if blast_files:
            cmd.extend(["-b", *(str(path) for path in blast_files)])
        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=output_dir,
        )
        return (
            result.returncode,
            result.stdout + result.stderr,
            output_dir / f"{output_prefix}.svg",
        )


@pytest.fixture(scope="session")
def gbdraw_runner() -> GbdrawRunner:
    return GbdrawRunner()


@pytest.fixture
def stub_typed_request_export(monkeypatch: pytest.MonkeyPatch) -> None:
    """Avoid writing renderer output in request-plumbing tests."""

    import gbdraw.api.request_render as request_render_module

    monkeypatch.setattr(
        request_render_module,
        "save_figure_to",
        lambda *_args, output_dir=None, output_prefix=None, **_kwargs: [
            str(Path(output_dir or ".") / f"{output_prefix}.svg")
        ],
    )
