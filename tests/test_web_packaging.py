from __future__ import annotations

import base64
from collections.abc import Callable
import gzip
import hashlib
import html
import importlib.util
import json
import re
import shutil
import socket
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from PIL import Image

from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    SUPPORTED_SESSION_VERSIONS,
)
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
README_PATH = REPO_ROOT / "README.md"
ABOUT_PATH = REPO_ROOT / "docs" / "ABOUT.md"
CITATION_PATH = REPO_ROOT / "CITATION.cff"

PREPRINT_TITLE = "gbdraw: a genome diagram generator for microbes and organelles"
PREPRINT_DOI = "10.64898/2026.04.07.716863"
BROWSER_WHEEL_FORBIDDEN_PREFIXES = (
    "gbdraw/web/assets/",
    "gbdraw/web/gallery/",
    "gbdraw/web/js/",
    "gbdraw/web/presets/",
    "gbdraw/web/tutorial-data/",
    "gbdraw/web/vendor/",
    "gbdraw/web/wasm/",
)
BROWSER_WHEEL_FORBIDDEN_FILES = {
    "gbdraw/web/index.html",
    "gbdraw/web/open-source-notices.html",
}
BROWSER_WHEEL_ALLOWED_FILES = {
    "gbdraw/web/js/services/standalone-interactivity-assets.js",
}
BROWSER_WHEEL_REQUIRED_RUNTIME_DATA = {
    "gbdraw/data/color_palettes.toml",
    "gbdraw/data/config.toml",
    "gbdraw/web/js/services/standalone-interactivity-assets.js",
}
GALLERY_SESSION_FILES = [
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_basic_circular.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "tobacco-chloroplast.gbdraw-session.json",
    "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "lambda_basic_linear.gbdraw-session.json",
]
GALLERY_MULTI_RECORD_LINEAR_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
}
GALLERY_EDITOR_STATE_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "WSSV_genome_comparison.gbdraw-session.json",
}
GALLERY_LOSAT_CACHE_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
}
GALLERY_LOSAT_DERIVED_CACHE_SESSION_FILES: set[str] = set()


def _load_verify_module():
    module_path = REPO_ROOT / "tools" / "verify_gui_offline.py"
    spec = importlib.util.spec_from_file_location("verify_gui_offline", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load verification module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_session_prefix(path: Path, length: int = 256) -> str:
    with path.open("rb") as session_file:
        is_gzip = session_file.read(2) == b"\x1f\x8b"
    if is_gzip:
        with gzip.open(path, mode="rt", encoding="utf-8") as session_handle:
            return session_handle.read(length)
    with path.open(encoding="utf-8") as session_handle:
        return session_handle.read(length)


def _load_prepare_browser_wheel_module():
    module_path = REPO_ROOT / "tools" / "prepare_browser_wheel.py"
    spec = importlib.util.spec_from_file_location("prepare_browser_wheel", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Could not load browser wheel preparation module from {module_path}"
        )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_prepare_cloudflare_pages_module():
    module_path = REPO_ROOT / "tools" / "prepare_cloudflare_pages.py"
    spec = importlib.util.spec_from_file_location(
        "prepare_cloudflare_pages", module_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Could not load Cloudflare packaging module from {module_path}"
        )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _can_bind_loopback() -> bool:
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind(("127.0.0.1", 0))
    except OSError:
        return False
    return True


def _run_prepare_browser_wheel(*args: str) -> None:
    subprocess.run(
        [sys.executable, "tools/prepare_browser_wheel.py", *args],
        cwd=REPO_ROOT,
        check=True,
    )


def ensure_prepared_browser_wheel():
    verify_module = _load_verify_module()
    try:
        browser_wheel_path = (
            verify_module.BUILD_SUPPORT.validate_browser_wheel_prepared()
        )
    except (FileNotFoundError, RuntimeError):
        _run_prepare_browser_wheel()
        browser_wheel_path = (
            verify_module.BUILD_SUPPORT.validate_browser_wheel_prepared()
        )
    return verify_module, browser_wheel_path


def _gallery_svg_metadata(svg_source: str) -> dict[str, object]:
    metadata_match = re.search(
        r'<metadata(?P<attributes>[^>]*id="gbdraw-interactive-feature-metadata"[^>]*)>'
        r"(?P<payload>.*?)</metadata>",
        svg_source,
        re.S,
    )
    assert metadata_match, "missing interactive feature metadata"
    payload = html.unescape(metadata_match.group("payload"))
    if 'data-encoding="gzip-base64"' in metadata_match.group("attributes"):
        payload = gzip.decompress(base64.b64decode(payload)).decode("utf-8")
    return json.loads(payload)


def _assert_white_gallery_thumbnail(path: Path) -> None:
    image = Image.open(path).convert("RGB")
    assert image.size == (640, 360)

    width, height = image.size
    corners = [
        image.getpixel((0, 0)),
        image.getpixel((width - 1, 0)),
        image.getpixel((0, height - 1)),
        image.getpixel((width - 1, height - 1)),
    ]
    assert all(min(pixel) >= 245 for pixel in corners)

    border_pixels = []
    for x in range(width):
        border_pixels.append(image.getpixel((x, 0)))
        border_pixels.append(image.getpixel((x, height - 1)))
    for y in range(height):
        border_pixels.append(image.getpixel((0, y)))
        border_pixels.append(image.getpixel((width - 1, y)))
    average_luminance = sum(sum(pixel) / 3 for pixel in border_pixels) / len(
        border_pixels
    )
    assert average_luminance >= 245


@pytest.mark.browser
def test_web_offline_assets_can_be_prepared_for_packaging() -> None:
    verify_module, expected_wheel_path = ensure_prepared_browser_wheel()
    expected_wheel_name = "gbdraw-0.14.0b0-py3-none-any.whl"
    assert verify_module._parse_wheel_name() == expected_wheel_name
    assert expected_wheel_path.name == expected_wheel_name
    verify_module.assert_browser_wheel_is_not_recursive(expected_wheel_path)
    verify_module.BUILD_SUPPORT.inspect_browser_wheel_payload(expected_wheel_path)
    verify_module._assert_packaged_assets()


@pytest.mark.browser
def test_browser_wheel_excludes_hosted_web_assets() -> None:
    verify_module, browser_wheel_path = ensure_prepared_browser_wheel()
    with zipfile.ZipFile(browser_wheel_path) as zf:
        names = set(zf.namelist())

    forbidden = sorted(
        name
        for name in names
        if name in BROWSER_WHEEL_FORBIDDEN_FILES
        or (
            name not in BROWSER_WHEEL_ALLOWED_FILES
            and any(
                name.startswith(prefix) for prefix in BROWSER_WHEEL_FORBIDDEN_PREFIXES
            )
        )
        or (name.startswith("gbdraw/web/gbdraw-") and name.endswith(".whl"))
    )
    assert forbidden == []
    assert BROWSER_WHEEL_REQUIRED_RUNTIME_DATA <= names
    assert (
        browser_wheel_path.stat().st_size
        <= verify_module.BUILD_SUPPORT.BROWSER_WHEEL_MAX_BYTES
    )


def test_local_web_package_data_excludes_gallery_assets() -> None:
    build_support = _load_verify_module().BUILD_SUPPORT
    package_data_patterns = build_support.get_package_data_patterns(
        include_browser_wheel=True
    )

    assert all("web/gallery" not in pattern for pattern in package_data_patterns)
    assert "web/js/services/*.js" in package_data_patterns
    assert "web/tutorial-data/*.json" in package_data_patterns
    assert "web/tutorial-data/*/*.gb" in package_data_patterns
    assert "web/tutorial-data/*/*.gbk" in package_data_patterns
    assert "web/tutorial-data/*/*.gff3" in package_data_patterns
    assert "web/tutorial-data/*/*.tsv" in package_data_patterns
    assert (WEB_ROOT / "js" / "services" / "canonical-comparisons.js").is_file()
    manifest_in = (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "gbdraw/web/gallery" not in manifest_in
    assert "recursive-include gbdraw/web/tutorial-data *" in manifest_in
    assert "include tools/build_lambda_gff3_fixture.py" in manifest_in
    assert "gbdraw/web/gallery/" in build_support._BROWSER_WHEEL_FORBIDDEN_PREFIXES
    assert (
        "gbdraw/web/tutorial-data/" in build_support._BROWSER_WHEEL_FORBIDDEN_PREFIXES
    )


def test_index_links_to_open_source_notices() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "./open-source-notices.html" in index_html
    assert "Open Source Notices" in index_html


@pytest.mark.gallery
def test_index_links_to_hosted_interactive_gallery() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "Interactive Gallery" in index_html
    assert "./gallery/" not in index_html

    match = re.search(
        r'<a href="https://gbdraw\.app/gallery/"(?P<attrs>[^>]*)>Interactive Gallery</a>',
        index_html,
    )
    assert match is not None
    assert 'target="_blank"' in match.group("attrs")
    assert 'rel="noopener noreferrer"' in match.group("attrs")


def test_public_web_html_entrypoints_are_not_gitignored() -> None:
    required_html_paths = [
        "gbdraw/web/index.html",
        "gbdraw/web/open-source-notices.html",
        "gbdraw/web/gallery/index.html",
    ]
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")

    for path in required_html_paths:
        assert f"!{path}" in gitignore

    if shutil.which("git") and (REPO_ROOT / ".git").exists():
        for path in required_html_paths:
            result = subprocess.run(
                ["git", "check-ignore", "-q", path],
                cwd=REPO_ROOT,
                check=False,
            )
            assert result.returncode == 1, (
                f"{path} must be commit-visible for hosted builds"
            )


@pytest.mark.gallery
def test_gallery_csp_allows_same_origin_tutorial_media() -> None:
    gallery_index = (WEB_ROOT / "gallery" / "index.html").read_text(encoding="utf-8")
    assert "media-src 'self';" in gallery_index
    assert "img-src 'self' data:;" in gallery_index


def test_open_source_notices_are_generated() -> None:
    subprocess.run(
        [sys.executable, "tools/generate_open_source_notices.py", "--check"],
        cwd=REPO_ROOT,
        check=True,
    )


def test_open_source_notices_generator_does_not_require_runtime_dependencies() -> None:
    subprocess.run(
        [sys.executable, "-S", "tools/generate_open_source_notices.py", "--check"],
        cwd=REPO_ROOT,
        check=True,
    )


def test_open_source_notices_include_current_project_version() -> None:
    from gbdraw._build_support import read_project_version

    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    assert read_project_version() in notices_html


def test_open_source_notices_omit_internal_generation_details() -> None:
    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    assert "Distribution Summary" not in notices_html
    assert "Inventory Sources" not in notices_html
    assert "Project version:" not in notices_html
    assert "tools/generate_open_source_notices.py" not in notices_html


def test_open_source_notices_include_local_pyodide_wheels() -> None:
    from gbdraw._web_assets import PYODIDE_LOCAL_WHEELS

    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    for wheel_path in PYODIDE_LOCAL_WHEELS:
        assert wheel_path.as_posix() in notices_html


def test_web_config_pyodide_local_wheels_match_shared_manifest() -> None:
    from gbdraw._web_assets import PYODIDE_LOCAL_WHEELS

    config_js = (WEB_ROOT / "js" / "config.js").read_text(encoding="utf-8")
    match = re.search(
        r"export const PYODIDE_LOCAL_WHEELS\s*=\s*\[(.*?)\];", config_js, re.DOTALL
    )
    assert match is not None
    configured_wheels = tuple(
        Path(raw.lstrip("./"))
        for raw in re.findall(r"""["']([^"']+\.whl)["']""", match.group(1))
    )
    assert configured_wheels == PYODIDE_LOCAL_WHEELS


def test_index_uses_title_logo_separately_from_icon_assets() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "./assets/gbdraw-logo-title.png" in index_html
    assert (
        '<link rel="icon" href="./assets/gbdraw-logo.svg" type="image/svg+xml">'
        in index_html
    )


def test_meta_csp_omits_frame_ancestors_header_only_directive() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    gallery_html = (GALLERY_ROOT / "index.html").read_text(encoding="utf-8")
    assert "frame-ancestors" not in index_html
    assert "frame-ancestors" not in notices_html
    assert "frame-ancestors" not in gallery_html


@pytest.mark.gallery
def test_interactive_gallery_shell_is_sandboxed() -> None:
    gallery_html = (GALLERY_ROOT / "index.html").read_text(encoding="utf-8")

    assert "default-src 'self';" in gallery_html
    assert "script-src 'self';" in gallery_html
    assert "style-src 'self';" in gallery_html
    assert "frame-src 'self';" in gallery_html
    assert '<script type="module" src="./gallery.js"></script>' in gallery_html
    assert 'sandbox="allow-scripts allow-downloads"' in gallery_html
    assert "allow-same-origin" not in gallery_html


@pytest.mark.gallery
def test_interactive_gallery_examples_are_wired() -> None:
    expected_ids = [
        "HmmtDNA_basic_circular",
        "lambda_basic_linear",
        "HmmtDNA_ATskew",
        "tobacco-chloroplast",
        "Vnig_TUMSAT-TG-2018",
        "hepatoplasmataceae_collinear",
        "vibrio-harveyi-group-collinear",
        "hepatoplasmataceae_orthogroup",
        "BGC0000708-BGC0000713",
        "majanivirus_orthogroup",
        "WSSV_genome_comparison",
    ]
    examples = json.loads((GALLERY_ROOT / "examples.json").read_text(encoding="utf-8"))
    expected_tags = {
        "HmmtDNA_basic_circular": ["Circular", "Interactive SVG"],
        "lambda_basic_linear": ["Linear", "Interactive SVG"],
        "HmmtDNA_ATskew": ["Circular", "Interactive SVG"],
        "tobacco-chloroplast": ["Circular", "Interactive SVG"],
        "Vnig_TUMSAT-TG-2018": ["Circular", "Multi-record", "Interactive SVG"],
        "hepatoplasmataceae_collinear": [
            "Linear",
            "Collinear groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "vibrio-harveyi-group-collinear": [
            "Linear",
            "Multi-record",
            "Collinear groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "hepatoplasmataceae_orthogroup": [
            "Linear",
            "Similarity groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "BGC0000708-BGC0000713": [
            "Linear",
            "Similarity groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "majanivirus_orthogroup": [
            "Linear",
            "Similarity groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "WSSV_genome_comparison": ["Circular", "LOSAT", "Interactive SVG"],
    }

    assert [entry["id"] for entry in examples] == expected_ids
    assert [entry["displayOrder"] for entry in examples] == sorted(
        entry["displayOrder"] for entry in examples
    )
    assert [entry["title"] for entry in examples] == [
        "Human mitochondrial genome: first circular figure",
        "Lambda phage: first linear figure",
        "Human mitochondrial genome (AT skew)",
        "<i>Nicotiana tabacum</i> chloroplast genome regions",
        "<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018",
        "Hepatoplasmataceae collinear protein-match blocks",
        "<i>Vibrio</i> Harveyi group multi-record collinearity",
        "Hepatoplasmataceae CDS protein-similarity links",
        "Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.",
        "Majanivirus CDS protein-similarity links",
        "White spot syndrome virus nucleotide-similarity rings",
    ]
    for entry in examples:
        assert entry["title"]
        assert entry["description"]
        assert "difficulty" not in entry
        assert entry["workflow"]
        assert entry["inputSummary"]
        assert "estimatedTime" not in entry
        assert isinstance(entry["displayOrder"], int)
        assert entry["commandKind"] in {"runnable", "provenance"}
        assert entry["commandNote"]
        assert not entry.get("interactiveStep")
        assert entry["tags"] == expected_tags[entry["id"]]
        assert entry["command"].startswith("gbdraw ")
        assert "interactive-svg" not in entry["command"]
        assert entry["fileSizeLabel"]
        assert entry["sourceNote"]
        assert entry["featureSources"]
        assert entry["svg"].startswith("./examples/")
        assert entry["session"].startswith("./sessions/")
        assert entry["thumbnail"].startswith("./thumbnails/")
        assert entry["tutorial"].startswith("./tutorials/")
        assert entry["tutorialStatus"] == "ready"
        assert entry["sourceSession"].startswith("gbdraw/web/gallery/sessions/")
        assert entry["sourceOutput"].startswith("gbdraw/web/gallery/examples/")
        assert entry["sourceFigure"].startswith("gbdraw/web/gallery/sources/")

        svg_path = GALLERY_ROOT / entry["svg"].removeprefix("./")
        session_path = GALLERY_ROOT / entry["session"].removeprefix("./")
        source_figure_path = REPO_ROOT / entry["sourceFigure"]
        thumbnail_path = GALLERY_ROOT / entry["thumbnail"].removeprefix("./")
        assert svg_path.exists()
        assert session_path.exists()
        assert source_figure_path.exists()
        assert thumbnail_path.exists()

        svg_source = svg_path.read_text(encoding="utf-8")
        session_prefix = _read_session_prefix(session_path)
        thumbnail_header = thumbnail_path.read_bytes()[:16]

        assert svg_path.stat().st_size > 1024
        assert session_path.stat().st_size > 1024
        assert session_path.stat().st_size < 100_000_000
        assert source_figure_path.stat().st_size > 1024
        assert thumbnail_path.stat().st_size > 1024
        assert '"format":"gbdraw-session"' in session_prefix
        version_match = re.search(r'"version":(\d+)', session_prefix)
        assert version_match is not None
        assert int(version_match.group(1)) in SUPPORTED_SESSION_VERSIONS
        assert "gbdraw-gallery-interactive-script" not in svg_source
        assert "data-gbdraw-gallery" not in svg_source
        assert "window.parent" not in svg_source
        assert "parent.postMessage" not in svg_source
        assert "window.top" not in svg_source
        assert "window.opener" not in svg_source

        if entry.get("svgType") == "static":
            assert 'data-gbdraw-interactive-svg="true"' not in svg_source
            assert "gbdraw-interactive-feature-script" not in svg_source
        else:
            assert 'data-gbdraw-interactive-svg="true"' in svg_source
            assert "gbdraw-interactive-feature-metadata" in svg_source
            assert "gbdraw-interactive-feature-script" in svg_source
            assert 'data-popup-mode="rich"' in svg_source
            assert "data-gbdraw-original-viewbox" in svg_source
            payload = _gallery_svg_metadata(svg_source)
            assert payload["schema"] == 3
            assert len(payload["items"]) == 1
            item = payload["items"][0]
            biological_features = item["biologicalFeatures"]
            rendered_features = item["features"]
            assert biological_features
            assert rendered_features
            assert any(feature.get("qualifiers") for feature in biological_features)
            assert all(
                isinstance(feature.get("start"), int)
                and isinstance(feature.get("end"), int)
                for feature in biological_features
            )
            assert all(
                "location_parts" not in feature
                or (
                    isinstance(feature["location_parts"], list)
                    and feature["location_parts"]
                    and all(
                        isinstance(part, dict) for part in feature["location_parts"]
                    )
                )
                for feature in biological_features
            )
            sequence_sources = item.get("sequenceSources", [])
            assert all(
                feature.get("nucleotide_sequence")
                or (
                    isinstance(feature.get("sequenceSourceIndex"), int)
                    and 0 <= feature["sequenceSourceIndex"] < len(sequence_sources)
                    and sequence_sources[feature["sequenceSourceIndex"]].get("sequence")
                )
                for feature in biological_features
            )
            biological_keys = {
                (feature["recordKey"], feature["biologicalFeatureId"])
                for feature in biological_features
            }
            assert all(
                (feature["recordKey"], feature["biologicalFeatureId"])
                in biological_keys
                for feature in rendered_features
            )
            for group in item["orthogroups"]:
                assert group["members"]
                assert all(
                    (member["recordKey"], member["biologicalFeatureId"])
                    in biological_keys
                    for member in group["members"]
                )
            for match in item["comparisonMatches"]:
                for role in ("query", "subject"):
                    record_key = match.get(f"{role}RecordKey")
                    feature_id = match.get(f"{role}BiologicalFeatureId")
                    if record_key or feature_id:
                        assert (record_key, feature_id) in biological_keys

        assert thumbnail_header.startswith(b"RIFF")
        assert b"WEBP" in thumbnail_header
        _assert_white_gallery_thumbnail(thumbnail_path)

    provenance = [entry for entry in examples if entry["commandKind"] == "provenance"]
    assert [entry["id"] for entry in provenance] == ["WSSV_genome_comparison"]
    assert "not directly runnable" in provenance[0]["commandNote"]
    collinear = next(
        entry for entry in examples if entry["id"] == "hepatoplasmataceae_collinear"
    )
    assert collinear["command"].count("--losatp_threads") == 1


@pytest.mark.gallery
def test_runnable_gallery_support_downloads_exist() -> None:
    expected_local_downloads = {
        "HmmtDNA_ATskew": {"./files/HmmtDNA_qualifier_priority.tsv"},
        "tobacco-chloroplast": {
            "./files/chloroplast_specific_table.tsv",
            "./files/qualifier_priority.tsv",
            "./files/nicotiana-tabacum-regions.tsv",
        },
        "BGC0000708-BGC0000713": {
            "./files/BGC0000708-BGC0000713_default_colors.tsv",
            "./files/BGC0000708-BGC0000713_specific_colors.tsv",
            "./files/BGC0000708-BGC0000713_qualifier_priority.tsv",
        },
        "majanivirus_orthogroup": {
            "./files/modified_default_colors.tsv",
            "./files/majani_custom_color_table.tsv",
        },
    }

    for example_id, hrefs in expected_local_downloads.items():
        tutorial = json.loads(
            (GALLERY_ROOT / "tutorials" / f"{example_id}.json").read_text(
                encoding="utf-8"
            )
        )
        tutorial_hrefs = {
            download.get("href") for download in tutorial.get("downloads", [])
        }
        assert hrefs <= tutorial_hrefs
        for href in hrefs:
            assert (GALLERY_ROOT / href.removeprefix("./")).is_file()


def test_index_includes_preprint_citation() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "How to cite" in index_html
    assert PREPRINT_TITLE in index_html
    assert PREPRINT_DOI in index_html


@pytest.mark.gallery
def test_gallery_sessions_ship_resumable_state_without_duplicate_files(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    observed_losat_cache_sessions = set()
    observed_losat_derived_cache_sessions = set()

    for session_name in GALLERY_SESSION_FILES:
        session_path = GALLERY_ROOT / "sessions" / session_name
        session = load_cached_gallery_session(session_path)
        results = session.get("results", [])
        editor_state = session.get("editorState", {})
        feature_catalog = editor_state.get("featureCatalog", {})
        catalog_items = feature_catalog.get("items", [])
        svg_text = "\n".join(result.get("content", "") for result in results)
        pairwise_ids = set(
            re.findall(r"data-gbdraw-pairwise-match-id=[\"']([^\"']+)[\"']", svg_text)
        )
        collinearity_ids = set(
            re.findall(r"data-collinearity-block-id=[\"']([^\"']+)[\"']", svg_text)
        )

        assert session.get("version") == CURRENT_SESSION_VERSION, session_name
        assert (
            session.get("renderRequest", {}).get("schema") == CANONICAL_REQUEST_SCHEMA
        ), session_name
        assert (
            session.get("proteinIdentityManifest", {}).get("schema")
            == PROTEIN_IDENTITY_MANIFEST_SCHEMA
        ), session_name
        assert "files" not in session, session_name
        assert results, session_name
        assert feature_catalog.get("schema") == 3, session_name
        assert [
            (item.get("resultIndex"), item.get("resultName")) for item in catalog_items
        ] == [
            (result_index, result.get("name"))
            for result_index, result in enumerate(results)
        ], session_name
        feature_state = session.get("features")
        if isinstance(feature_state, dict):
            assert "extractedFeatures" not in feature_state, session_name
        orthogroup_state = session.get("orthogroupState")
        if isinstance(orthogroup_state, dict):
            assert "groups" not in orthogroup_state, session_name
        assert "orthogroupState" in session, session_name

        for result_index, result in enumerate(results):
            assert (
                not str(result.get("name") or "").lower().endswith(".interactive.svg")
            ), session_name
            content = result.get("content", "")
            assert "<svg" in content, session_name
            assert "gbdraw-interactive-feature-metadata" not in content, session_name

            item = catalog_items[result_index]
            biological_features = item.get("biologicalFeatures", [])
            rendered_features = item.get("features", [])
            biological_keys = {
                (
                    feature.get("recordKey"),
                    feature.get("biologicalFeatureId"),
                )
                for feature in biological_features
            }
            assert biological_keys, session_name
            assert rendered_features, session_name
            assert all(
                (
                    feature.get("recordKey"),
                    feature.get("biologicalFeatureId"),
                )
                in biological_keys
                for feature in rendered_features
            ), session_name
            catalog_svg_ids = {
                feature.get("svgId")
                for feature in rendered_features
                if feature.get("svgId")
            }
            result_svg_ids = set(
                re.findall(
                    r"data-gbdraw-feature-id=[\"']([^\"']+)[\"']",
                    content,
                )
            )
            assert result_svg_ids, session_name
            assert result_svg_ids <= catalog_svg_ids, session_name

        seen_resource_payloads: dict[tuple[int, bytes], str] = {}
        for resource_id, resource in session.get("resources", {}).items():
            assert resource.get("encoding") == "base64", session_name
            payload = base64.b64decode(resource.get("data", ""), validate=True)
            payload_identity = (
                len(payload),
                hashlib.sha256(payload).digest(),
            )
            assert payload_identity not in seen_resource_payloads, (
                session_name,
                resource_id,
                seen_resource_payloads.get(payload_identity),
            )
            seen_resource_payloads[payload_identity] = resource_id

        if session_name in GALLERY_EDITOR_STATE_SESSION_FILES:
            assert session.get("editorState"), session_name
        losat_cache_entries = session.get("losatCache", {}).get("entries", [])
        if losat_cache_entries:
            observed_losat_cache_sessions.add(session_name)
            for entry in losat_cache_entries:
                expected_schema = (
                    PROTEIN_LOSAT_CACHE_SCHEMA
                    if entry.get("identityKind") == "protein"
                    else NUCLEOTIDE_LOSAT_CACHE_SCHEMA
                )
                assert entry.get("schema") == expected_schema, session_name

        losat_derived_cache_entries = session.get("losatDerivedCache", {}).get(
            "entries", []
        )
        if losat_derived_cache_entries:
            observed_losat_derived_cache_sessions.add(session_name)
            assert all(
                entry.get("schema") == LOSAT_DERIVED_CACHE_SCHEMA
                for entry in losat_derived_cache_entries
            ), session_name
        if session_name in GALLERY_MULTI_RECORD_LINEAR_SESSION_FILES:
            assert pairwise_ids or collinearity_ids, session_name

    assert observed_losat_cache_sessions == GALLERY_LOSAT_CACHE_SESSION_FILES
    assert (
        observed_losat_derived_cache_sessions
        == GALLERY_LOSAT_DERIVED_CACHE_SESSION_FILES
    )


@pytest.mark.gallery
def test_tobacco_gallery_session_keeps_chloroplast_region_annotations(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    session = load_cached_gallery_session(
        GALLERY_ROOT / "sessions" / "tobacco-chloroplast.gbdraw-session.json"
    )
    annotation_sets = (
        session.get("renderRequest", {})
        .get("diagramOptions", {})
        .get("annotations", {})
        .get("sets", [])
    )
    assert len(annotation_sets) == 1
    assert annotation_sets[0]["id"] == "plastome_regions"
    annotations = annotation_sets[0]["annotations"]
    assert [annotation["label"] for annotation in annotations] == [
        "LSC",
        "IRb",
        "SSC",
        "IRa",
    ]
    actual_spans = [
        (annotation["target"]["start"], annotation["target"]["end"])
        for annotation in annotations
    ]
    assert actual_spans == [
        (1, 86686),
        (86687, 112029),
        (112030, 130600),
        (130601, 155943),
    ]


@pytest.mark.gallery
def test_vnig_gallery_session_multirecord_positions_are_restoreable(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    session = load_cached_gallery_session(
        GALLERY_ROOT / "sessions" / "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz"
    )
    expected_positions = ["#1@1", "#2@1", "#3@2", "#4@2", "#5@2", "#6@2"]
    config_positions = (
        session.get("config", {}).get("adv", {}).get("multi_record_positions")
    )
    if isinstance(config_positions, list) and config_positions:
        actual_positions = [
            f"{entry.get('selector')}@{entry.get('row')}" for entry in config_positions
        ]
    else:
        args = session.get("cliInvocation", {}).get("args", [])
        actual_positions = [
            str(args[index + 1])
            for index, arg in enumerate(args[:-1])
            if arg == "--multi_record_position"
        ]

    assert actual_positions == expected_positions


def test_prepare_browser_wheel_refreshes_open_source_notices(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prepare_module = _load_prepare_browser_wheel_module()
    repo_root = tmp_path / "repo"
    web_root = repo_root / "gbdraw" / "web"
    web_root.mkdir(parents=True)
    expected_name = "gbdraw-0.14.0b0-py3-none-any.whl"
    calls: list[object] = []

    def fake_run(
        args: list[str], *, cwd: Path, env: dict[str, str], check: bool
    ) -> None:
        calls.append("wheel")
        assert cwd == repo_root
        assert check is True
        assert env["GBDRAW_BUILDING_BROWSER_WHEEL"] == "1"
        assert args[:4] == [sys.executable, "-m", "pip", "wheel"]
        assert "--no-deps" in args
        assert "--no-cache-dir" in args
        assert "--no-build-isolation" not in args
        assert args[-1] == str(repo_root)
        outdir = Path(args[args.index("--wheel-dir") + 1])
        outdir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(outdir / expected_name, "w") as zf:
            zf.writestr("gbdraw/__init__.py", "")

    def fake_validate_browser_wheel_prepared() -> Path:
        calls.append("validate")
        return web_root / expected_name

    build_support = SimpleNamespace(
        BROWSER_WHEEL_BUILD_ENV="GBDRAW_BUILDING_BROWSER_WHEEL",
        expected_browser_wheel_name=lambda: expected_name,
        refresh_open_source_notices=lambda: calls.append("notices"),
        generate_cache_bust_token=lambda: "cache-token",
        update_browser_wheel_config=lambda **kwargs: calls.append(("config", kwargs)),
        validate_browser_wheel_prepared=fake_validate_browser_wheel_prepared,
    )
    monkeypatch.setattr(prepare_module, "REPO_ROOT", repo_root)
    monkeypatch.setattr(prepare_module, "WEB_ROOT", web_root)
    monkeypatch.setattr(
        prepare_module, "_load_build_support_module", lambda: build_support
    )
    monkeypatch.setattr(prepare_module.subprocess, "run", fake_run)

    assert prepare_module.prepare_browser_wheel(refresh_cache_bust=True) == 0
    assert (web_root / expected_name).exists()
    assert calls == [
        "notices",
        "wheel",
        ("config", {"wheel_name": expected_name, "cache_bust": "cache-token"}),
        "validate",
    ]


@pytest.mark.browser
def test_cloudflare_bundle_includes_google_analytics_and_hosted_notice(
    tmp_path: Path,
) -> None:
    from gbdraw._build_support import read_project_version

    verify_module, _ = ensure_prepared_browser_wheel()
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    output_root = tmp_path / "cloudflare-pages"
    commit_sha = "abcdef1234567890abcdef1234567890abcdef12"
    remote_base = (
        "https://raw.githubusercontent.com/satoshikawato/gbdraw/"
        f"{commit_sha}/gbdraw/web/"
    )
    bundle_path = cloudflare_module.build_cloudflare_pages_bundle(
        output_root=output_root,
        commit_sha=commit_sha,
    )

    missing_tutorial_assets = [
        path.as_posix()
        for path in verify_module.REQUIRED_TUTORIAL_DATA_FILES
        if not (bundle_path / path).is_file()
    ]
    assert missing_tutorial_assets == []

    index_html = (bundle_path / "index.html").read_text(encoding="utf-8")
    assert "https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y" in index_html
    assert "gtag('config', 'G-GG6JMKM02Y');" in index_html
    assert "static.cloudflareinsights.com" not in index_html
    assert "cloudflareinsights.com" not in index_html
    assert "Hosted Site Analytics" in index_html
    assert "uses Google Analytics 4 for aggregate page-usage metrics" in index_html
    assert (
        "Uploaded genome files and generated diagrams are still processed locally in your browser"
        in index_html
    )
    assert (
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://*.googletagmanager.com;"
        in index_html
    )
    assert (
        "img-src 'self' data: blob: https://*.google-analytics.com "
        "https://*.googletagmanager.com;"
    ) in index_html
    assert (
        "connect-src 'self' https://*.google-analytics.com "
        "https://*.analytics.google.com https://*.googletagmanager.com;"
    ) in index_html
    assert "GOOGLE_ANALYTICS_SCRIPT" not in index_html
    assert "GOOGLE_ANALYTICS_NOTICE" not in index_html
    assert "GBDRAW_HOSTED_BUILD_LABEL" not in index_html
    assert f"Version: v{read_project_version()}+abcdef1" in index_html
    assert f'title="Commit {commit_sha}"' in index_html
    headers = (bundle_path / "_headers").read_text(encoding="utf-8")
    assert "Cross-Origin-Opener-Policy: same-origin" in headers
    assert "Cross-Origin-Embedder-Policy: require-corp" in headers
    assert "Cross-Origin-Resource-Policy: same-origin" in headers
    assert "Content-Security-Policy: frame-ancestors 'none'" in headers
    assert "/gallery/examples/*" in headers
    assert "! Content-Security-Policy" in headers
    assert "frame-ancestors 'self'" in headers
    assert "gallery/media/**/*" in cloudflare_module.GALLERY_REMOTE_ASSET_PATTERNS
    assert (
        "gallery/sessions/*.gbdraw-session.json.gz"
        in cloudflare_module.GALLERY_REMOTE_ASSET_PATTERNS
    )
    remote_assets = json.loads(
        (bundle_path / "gallery" / "remote-assets.json").read_text(encoding="utf-8")
    )
    assert (
        remote_assets["gallery/examples/Vnig_TUMSAT-TG-2018.svg"]
        == f"{remote_base}gallery/examples/Vnig_TUMSAT-TG-2018.svg"
    )
    assert (
        "gallery/sessions/Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz"
        not in remote_assets
    )
    assert (
        remote_assets[
            "gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
        ]
        == f"{remote_base}gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
    )
    assert (
        remote_assets["gallery/examples/vibrio-harveyi-group-collinear.svg"]
        == f"{remote_base}gallery/examples/vibrio-harveyi-group-collinear.svg"
    )
    assert all("/main/" not in url for url in remote_assets.values())
    assert not (
        bundle_path / "gallery" / "examples" / "Vnig_TUMSAT-TG-2018.svg"
    ).exists()
    assert not (
        bundle_path / "gallery" / "examples" / "vibrio-harveyi-group-collinear.svg"
    ).exists()
    assert (
        bundle_path
        / "gallery"
        / "sessions"
        / "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz"
    ).exists()
    assert not (
        bundle_path
        / "gallery"
        / "sessions"
        / "vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
    ).exists()
    assert (
        bundle_path / "gallery" / "examples" / "majanivirus_orthogroup.svg"
    ).exists()


def test_cloudflare_gallery_remote_base_rejects_mutable_refs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    monkeypatch.delenv(cloudflare_module.GALLERY_REMOTE_BASE_ENV, raising=False)
    monkeypatch.setenv(cloudflare_module.GALLERY_REMOTE_REF_ENV, "main")

    with pytest.raises(RuntimeError, match="full 40-character commit SHA"):
        cloudflare_module._default_gallery_remote_base()


def test_cloudflare_gallery_remote_base_uses_resolved_commit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    commit_sha = "1234567890abcdef1234567890abcdef12345678"
    monkeypatch.delenv(cloudflare_module.GALLERY_REMOTE_BASE_ENV, raising=False)
    monkeypatch.delenv(cloudflare_module.GALLERY_REMOTE_REF_ENV, raising=False)
    monkeypatch.setattr(
        cloudflare_module,
        "_load_stamp_web_build_module",
        lambda: SimpleNamespace(resolve_commit_sha=lambda: commit_sha),
    )

    assert cloudflare_module._default_gallery_remote_base() == (
        "https://raw.githubusercontent.com/satoshikawato/gbdraw/"
        f"{commit_sha}/gbdraw/web/"
    )


def test_cloudflare_prepare_refreshes_gallery_only_when_requested(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    calls: list[tuple[str, object]] = []
    output_root = tmp_path / "cloudflare-pages"

    monkeypatch.setattr(
        cloudflare_module,
        "_load_prepare_browser_wheel_module",
        lambda: SimpleNamespace(
            prepare_browser_wheel=lambda refresh_cache_bust=False: calls.append(
                ("wheel", refresh_cache_bust)
            )
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_build_support_module",
        lambda: SimpleNamespace(
            refresh_open_source_notices=lambda: calls.append(("notices", None))
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_refresh_gallery_sessions_module",
        lambda: pytest.fail("Gallery refresh should be opt-in for Cloudflare Pages."),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "build_cloudflare_pages_bundle",
        lambda *, output_root=cloudflare_module.DEFAULT_OUTPUT_ROOT, google_analytics_measurement_id=cloudflare_module.DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID, gallery_remote_base=None: (
            calls.append(("bundle", output_root)) or output_root
        ),
    )

    assert (
        cloudflare_module.prepare_cloudflare_pages(output_root=output_root)
        == output_root
    )
    assert calls == [("wheel", False), ("notices", None), ("bundle", output_root)]

    calls.clear()
    gallery_module = SimpleNamespace(
        refresh_gallery_sessions=lambda: calls.append(("gallery", "sessions")),
        prepare_gallery_assets=lambda: calls.append(("gallery", "assets")),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_refresh_gallery_sessions_module",
        lambda: gallery_module,
    )

    assert (
        cloudflare_module.prepare_cloudflare_pages(
            refresh_cache_bust=True,
            refresh_gallery_sessions=True,
            output_root=output_root,
        )
        == output_root
    )
    assert calls == [
        ("gallery", "sessions"),
        ("gallery", "assets"),
        ("wheel", True),
        ("notices", None),
        ("bundle", output_root),
    ]


def test_cloudflare_prepare_refreshes_open_source_notices_before_copy(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    calls: list[str] = []
    output_root = tmp_path / "cloudflare-pages"

    monkeypatch.setattr(
        cloudflare_module,
        "_load_prepare_browser_wheel_module",
        lambda: SimpleNamespace(
            prepare_browser_wheel=lambda refresh_cache_bust=False: calls.append("wheel")
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_build_support_module",
        lambda: SimpleNamespace(
            refresh_open_source_notices=lambda: calls.append("notices")
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "build_cloudflare_pages_bundle",
        lambda *, output_root=cloudflare_module.DEFAULT_OUTPUT_ROOT, google_analytics_measurement_id=cloudflare_module.DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID, gallery_remote_base=None: (
            calls.append("copy") or output_root
        ),
    )

    assert (
        cloudflare_module.prepare_cloudflare_pages(output_root=output_root)
        == output_root
    )
    assert calls == [
        "wheel",
        "notices",
        "copy",
    ]


def test_wrangler_uses_cloudflare_bundle_directory() -> None:
    wrangler_toml = (REPO_ROOT / "wrangler.toml").read_text(encoding="utf-8")
    assert 'main = "./gbdraw/web/cloudflare-worker.js"' in wrangler_toml
    assert 'directory = "./dist/cloudflare-pages"' in wrangler_toml
    assert 'binding = "ASSETS"' in wrangler_toml
    assert 'not_found_handling = "single-page-application"' in wrangler_toml
    assert '"/gallery/*"' in wrangler_toml


def test_project_docs_and_citation_metadata_include_preprint_doi() -> None:
    readme = README_PATH.read_text(encoding="utf-8")
    assert PREPRINT_DOI in readme
    assert "./gbdraw/web/assets/gbdraw-logo-title.png" in readme
    assert PREPRINT_DOI in ABOUT_PATH.read_text(encoding="utf-8")
    citation_cff = CITATION_PATH.read_text(encoding="utf-8")
    assert PREPRINT_DOI in citation_cff
    assert "preferred-citation:" in citation_cff


def test_conda_build_prepares_browser_wheel_before_install() -> None:
    build_sh = (REPO_ROOT / "recipe" / "build.sh").read_text(encoding="utf-8")
    meta_yaml = (REPO_ROOT / "recipe" / "meta.yaml").read_text(encoding="utf-8")

    prepare_index = build_sh.index(
        "$PYTHON tools/prepare_browser_wheel.py --no-build-isolation"
    )
    install_index = build_sh.index(
        "$PYTHON -m pip install . --no-deps --ignore-installed -vv"
    )
    assert prepare_index < install_index
    assert "python-build" not in meta_yaml
    assert re.search(r"^\s+- setuptools\s*$", meta_yaml, re.MULTILINE)
    assert re.search(r"^\s+- wheel\s*$", meta_yaml, re.MULTILINE)


def test_hosted_web_build_refreshes_gallery_sessions_before_copy() -> None:
    deploy_yml = (REPO_ROOT / ".github" / "workflows" / "deploy_web.yml").read_text(
        encoding="utf-8"
    )
    cloudflare_source = (REPO_ROOT / "tools" / "prepare_cloudflare_pages.py").read_text(
        encoding="utf-8"
    )

    assert 'python -m pip install -e ".[dev]"' in deploy_yml
    refresh_index = deploy_yml.index("python tools/refresh_gallery_sessions.py")
    copy_index = deploy_yml.index("cp -r gbdraw/web/* public/")
    stamp_index = deploy_yml.index("python tools/stamp_web_build.py public")
    assert refresh_index < copy_index
    assert copy_index < stamp_index
    assert "refresh_gallery_sessions: bool = False" in cloudflare_source
    assert '"--refresh-gallery"' in cloudflare_source
    assert "refresh_gallery_sessions=args.refresh_gallery" in cloudflare_source
    assert (
        "refresh_gallery_sessions_module.refresh_gallery_sessions()"
        in cloudflare_source
    )
    assert (
        "refresh_gallery_sessions_module.prepare_gallery_assets()" in cloudflare_source
    )


def test_conda_recipe_does_not_copy_entire_web_tree() -> None:
    build_sh = (REPO_ROOT / "recipe" / "build.sh").read_text(encoding="utf-8")

    assert "cp -r gbdraw/web/*" not in build_sh
    assert "gbdraw/web/gallery" not in build_sh
    assert "gbdraw/web/index.html" in build_sh
    assert "gbdraw/web/open-source-notices.html" in build_sh
    assert "gbdraw/web/gbdraw-*.whl" in build_sh
    assert (
        "for web_asset_dir in assets js presets tutorial-data vendor wasm" in build_sh
    )


def test_setup_commands_refresh_open_source_notices() -> None:
    setup_source = (REPO_ROOT / "setup.py").read_text(encoding="utf-8")

    assert "class build_py(_build_py):" in setup_source
    assert "class sdist(_sdist):" in setup_source
    assert setup_source.count("_build_support.refresh_open_source_notices()") == 2
    assert 'cmdclass={"build_py": build_py, "sdist": sdist}' in setup_source


@pytest.mark.slow
def test_build_py_copies_offline_gui_assets(tmp_path: Path) -> None:
    verify_module, _ = ensure_prepared_browser_wheel()
    build_root = tmp_path / "build_lib"
    subprocess.run(
        [sys.executable, "setup.py", "build_py", "--build-lib", str(build_root)],
        cwd=REPO_ROOT,
        check=True,
    )

    required = [
        build_root / "gbdraw" / "web" / "index.html",
        build_root / "gbdraw" / "web" / "open-source-notices.html",
        build_root / "gbdraw" / "web" / "assets" / "favicon.ico",
        build_root / "gbdraw" / "web" / "assets" / "gbdraw-logo.svg",
        build_root / "gbdraw" / "web" / "assets" / "gbdraw-logo-title.svg",
        build_root / "gbdraw" / "web" / "assets" / "gbdraw-logo-title.png",
        build_root / "gbdraw" / "web" / verify_module._parse_wheel_name(),
        build_root / "gbdraw" / "web" / "vendor" / "vue" / "vue.global.js",
        build_root
        / "gbdraw"
        / "web"
        / "vendor"
        / "tailwindcss"
        / "tailwindcss-play.js",
        build_root
        / "gbdraw"
        / "web"
        / "vendor"
        / "pyodide"
        / "v0.29.0"
        / "full"
        / "pyodide.js",
        build_root
        / "gbdraw"
        / "web"
        / "vendor"
        / "pyodide"
        / "v0.29.0"
        / "full"
        / "pyodide.asm.wasm",
        build_root
        / "gbdraw"
        / "web"
        / "vendor"
        / "browser_wasi_shim"
        / "dist"
        / "index.js",
        build_root
        / "gbdraw"
        / "web"
        / "vendor"
        / "phosphor-icons"
        / "regular"
        / "style.css",
        build_root / "gbdraw" / "web" / "js" / "workers" / "losat-threaded-worker.js",
        build_root
        / "gbdraw"
        / "web"
        / "js"
        / "workers"
        / "losat-wasi-thread-worker.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "record-discovery.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "record-options.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "linear-record-selector.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "right-drawer.js",
        build_root
        / "gbdraw"
        / "web"
        / "js"
        / "app"
        / "annotations"
        / "record-catalog.js",
        build_root
        / "gbdraw"
        / "web"
        / "js"
        / "app"
        / "annotations"
        / "record-selector.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "annotations" / "validation.js",
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat.wasm",
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat-threaded.wasm",
        *(
            build_root / "gbdraw" / "web" / path
            for path in verify_module.REQUIRED_TUTORIAL_DATA_FILES
        ),
        *(
            build_root / "gbdraw" / "web" / path
            for path in verify_module.REQUIRED_UI_FONT_FILES
        ),
        *(
            build_root / "gbdraw" / "web" / path
            for path in verify_module._parse_local_wheel_paths()
        ),
    ]
    missing = [
        str(path.relative_to(build_root)) for path in required if not path.exists()
    ]
    assert not missing, (
        "build_py did not copy required offline GUI assets:\n" + "\n".join(missing)
    )
    assert not (build_root / "gbdraw" / "web" / "gallery").exists()
    copied_wheels = sorted(
        path.name for path in (build_root / "gbdraw" / "web").glob("gbdraw-*.whl")
    )
    assert copied_wheels == [verify_module._parse_wheel_name()]


@pytest.mark.slow
def test_built_wheel_contains_offline_gui_assets(tmp_path: Path) -> None:
    if importlib.util.find_spec("build") is None:
        pytest.skip("python -m build is not available in this environment")
    if importlib.util.find_spec("wheel") is None:
        pytest.skip("wheel is not available in this environment")

    verify_module, _ = ensure_prepared_browser_wheel()
    dist_dir = tmp_path / "dist"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--no-isolation",
            "--outdir",
            str(dist_dir),
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    wheel_path = next(dist_dir.glob("gbdraw-*.whl"))
    assert wheel_path.name == "gbdraw-0.14.0b0-py3-none-any.whl"
    subprocess.run(
        [
            sys.executable,
            "tools/verify_gui_offline.py",
            "inspect-wheel",
            str(wheel_path),
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    verify_module.assert_embedded_browser_wheel_is_not_recursive(wheel_path)

    with zipfile.ZipFile(wheel_path) as outer_wheel:
        outer_names = outer_wheel.namelist()
        browser_wheel_member = f"gbdraw/web/{verify_module._parse_wheel_name()}"
        browser_wheels = sorted(
            name
            for name in outer_names
            if name.startswith("gbdraw/web/gbdraw-") and name.endswith(".whl")
        )
        gallery_members = sorted(
            name for name in outer_names if name.startswith("gbdraw/web/gallery/")
        )
        assert browser_wheels == [browser_wheel_member]
        assert gallery_members == []
        assert "gbdraw/web/js/app/record-discovery.js" in outer_names
        assert "gbdraw/web/js/app/record-options.js" in outer_names
        assert "gbdraw/web/js/app/linear-record-selector.js" in outer_names
        assert "gbdraw/web/js/app/right-drawer.js" in outer_names
        assert "gbdraw/web/js/services/canonical-comparisons.js" in outer_names
        assert "gbdraw/web/js/app/annotations/record-catalog.js" in outer_names
        assert "gbdraw/web/js/app/annotations/record-selector.js" in outer_names
        assert "gbdraw/web/js/app/annotations/validation.js" in outer_names
        assert {
            f"gbdraw/web/{path.as_posix()}"
            for path in verify_module.REQUIRED_TUTORIAL_DATA_FILES
        } <= set(outer_names)


@pytest.mark.slow
def test_built_sdist_contains_tutorial_data(tmp_path: Path) -> None:
    if importlib.util.find_spec("build") is None:
        pytest.skip("python -m build is not available in this environment")

    verify_module, _ = ensure_prepared_browser_wheel()
    dist_dir = tmp_path / "dist"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--sdist",
            "--no-isolation",
            "--outdir",
            str(dist_dir),
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    sdist_path = next(dist_dir.glob("gbdraw-*.tar.gz"))
    with tarfile.open(sdist_path, "r:gz") as sdist:
        names = set(sdist.getnames())
    for path in verify_module.REQUIRED_TUTORIAL_DATA_FILES:
        suffix = f"/gbdraw/web/{path.as_posix()}"
        assert any(name.endswith(suffix) for name in names), suffix
    assert any(name.endswith("/tools/build_lambda_gff3_fixture.py") for name in names)


def _run_offline_gui_browser_contract(contract: str) -> None:
    if importlib.util.find_spec("playwright") is None:
        pytest.skip("playwright is not available in this environment")
    if not _can_bind_loopback():
        pytest.skip("loopback sockets are not permitted in this environment")

    ensure_prepared_browser_wheel()
    subprocess.run(
        [
            sys.executable,
            "tools/verify_gui_offline.py",
            "smoke-test",
            "--contract",
            contract,
        ],
        cwd=REPO_ROOT,
        check=True,
    )


@pytest.mark.slow
@pytest.mark.browser
def test_offline_gui_initializes_without_external_network_requests() -> None:
    _run_offline_gui_browser_contract("offline-initialization")


@pytest.mark.slow
@pytest.mark.browser
def test_offline_gui_generate_uses_one_lazy_diagram_worker() -> None:
    _run_offline_gui_browser_contract("generate")


@pytest.mark.slow
@pytest.mark.browser
def test_offline_gui_helper_then_render_reuses_one_diagram_worker() -> None:
    _run_offline_gui_browser_contract("helper-before-render")


@pytest.mark.slow
@pytest.mark.browser
def test_offline_gui_palette_preview_behavior() -> None:
    _run_offline_gui_browser_contract("palette-preview")


@pytest.mark.slow
@pytest.mark.browser
def test_offline_gui_exports_svg_png_pdf_and_interactive_svg() -> None:
    _run_offline_gui_browser_contract("exports")


@pytest.mark.slow
@pytest.mark.browser
def test_offline_gui_linear_losat_generation_populates_cache() -> None:
    _run_offline_gui_browser_contract("linear-losat")
