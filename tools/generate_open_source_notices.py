#!/usr/bin/env python3
from __future__ import annotations

import argparse
import difflib
import html
import sys
import zipfile
from dataclasses import dataclass
from email import policy
from email.message import Message
from email.parser import Parser
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw._build_support import read_project_version  # noqa: E402
from gbdraw._web_assets import (  # noqa: E402
    PYODIDE_BUNDLED_WHEELS,
    WEB_ASSET_NOTICES,
    WEB_ROOT,
)


NOTICES_PATH = WEB_ROOT / "open-source-notices.html"

WHEEL_LICENSE_OVERRIDES = {
    "bcbio-gff": "Biopython License",
    "biopython": "Biopython License",
    "fonttools": "MIT",
    "micropip": "MPL-2.0",
    "numpy": "BSD-3-Clause",
    "pandas": "BSD-3-Clause",
    "python-dateutil": "Apache-2.0 OR BSD-3-Clause",
    "pytz": "MIT",
    "six": "MIT",
    "svgwrite": "MIT",
    "tzdata": "Apache-2.0",
}


def _read_text_file(path: Path) -> str:
    return path.read_text(encoding="utf-8").strip()


def _read_wheel_member_text(relative_wheel_path: Path, member_suffix: str) -> str:
    wheel_path = WEB_ROOT / relative_wheel_path
    with zipfile.ZipFile(wheel_path) as zf:
        matches = sorted(name for name in zf.namelist() if name.endswith(member_suffix))
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one license member ending with {member_suffix!r} in "
                f"{wheel_path.relative_to(REPO_ROOT)}, found {len(matches)}"
            )
        return zf.read(matches[0]).decode("utf-8", "replace").strip()


LICENSE_TEXTS = {
    "license-mit": (
        "MIT License",
        lambda: _read_text_file(REPO_ROOT / "LICENSE.txt"),
    ),
    "license-apache-2.0": (
        "Apache License 2.0",
        lambda: _read_wheel_member_text(
            Path("vendor") / "pyodide" / "v0.29.0" / "full" / "tzdata-2025.2-py2.py3-none-any.whl",
            "licenses/LICENSE_APACHE",
        ),
    ),
    "license-mpl-2.0": (
        "Mozilla Public License 2.0",
        lambda: _read_wheel_member_text(
            Path("vendor") / "pyodide" / "v0.29.0" / "full" / "micropip-0.11.0-py3-none-any.whl",
            ".dist-info/licenses/LICENSE",
        ),
    ),
    "license-ofl-1.1": (
        "SIL Open Font License 1.1",
        lambda: _read_text_file(REPO_ROOT / "LICENSE_LIBERATION_FONTS.txt"),
    ),
    "license-bsd-3-clause": (
        "BSD 3-Clause License",
        lambda: _read_wheel_member_text(
            Path("vendor")
            / "pyodide-wheels"
            / "pandas-2.3.2-cp313-cp313-pyodide_2025_0_wasm32.whl",
            ".dist-info/LICENSE",
        ),
    ),
}

LICENSE_ANCHOR_LABELS = {
    "license-mit": "MIT",
    "license-apache-2.0": "Apache-2.0",
    "license-mpl-2.0": "MPL-2.0",
    "license-ofl-1.1": "OFL-1.1",
    "license-bsd-3-clause": "BSD-3-Clause",
}


@dataclass(frozen=True)
class ComponentNotice:
    display_name: str
    version: str
    license_expression: str
    source_url: str
    bundled_path: str
    notice: str
    license_text_anchors: tuple[str, ...]


@dataclass(frozen=True)
class WheelNotice:
    filename: str
    package_name: str
    version: str
    license_expression: str
    source_url: str
    bundled_path: str


def _escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def _external_link(url: str) -> str:
    if not url:
        return ""
    escaped = _escape(url)
    return f'<a href="{escaped}" target="_blank" rel="noopener noreferrer">{escaped}</a>'


def _slug(value: str) -> str:
    chars: list[str] = []
    previous_dash = False
    for char in value.lower():
        if char.isalnum():
            chars.append(char)
            previous_dash = False
        elif not previous_dash:
            chars.append("-")
            previous_dash = True
    return "".join(chars).strip("-") or "component"


def _license_anchors_for_expression(expression: str) -> tuple[str, ...]:
    normalized = expression.lower()
    anchors: list[str] = []
    if "mit" in normalized:
        anchors.append("license-mit")
    if "apache" in normalized:
        anchors.append("license-apache-2.0")
    if "mpl" in normalized or "mozilla public license" in normalized:
        anchors.append("license-mpl-2.0")
    if "ofl" in normalized or "open font license" in normalized:
        anchors.append("license-ofl-1.1")
    if "bsd" in normalized:
        anchors.append("license-bsd-3-clause")
    return tuple(dict.fromkeys(anchors))


def _license_html(expression: str, anchors: tuple[str, ...]) -> str:
    if not anchors:
        return _escape(expression)
    links = [
        f'<a href="#{_escape(anchor)}">{_escape(LICENSE_ANCHOR_LABELS.get(anchor, anchor))}</a>'
        for anchor in anchors
    ]
    return f"{_escape(expression)}. See {', '.join(links)}."


def _normalize_notice_text(text: str) -> str:
    return "\n".join(line.rstrip() for line in text.splitlines()).strip()


def _component_notices(project_version: str) -> tuple[ComponentNotice, ...]:
    project = ComponentNotice(
        display_name="gbdraw",
        version=project_version,
        license_expression="MIT",
        source_url="https://github.com/satoshikawato/gbdraw",
        bundled_path="gbdraw/web/",
        notice="Copyright (c) 2023 Satoshi Kawato.",
        license_text_anchors=("license-mit",),
    )
    assets = tuple(
        ComponentNotice(
            display_name=notice.display_name,
            version=notice.version,
            license_expression=notice.license_expression,
            source_url=notice.source_url,
            bundled_path=notice.bundled_path,
            notice=notice.notice,
            license_text_anchors=notice.license_text_anchors,
        )
        for notice in WEB_ASSET_NOTICES
    )
    return (project, *assets)


def _read_wheel_metadata(wheel_path: Path) -> Message:
    with zipfile.ZipFile(wheel_path) as zf:
        metadata_paths = sorted(name for name in zf.namelist() if name.endswith(".dist-info/METADATA"))
        if len(metadata_paths) != 1:
            raise RuntimeError(
                f"Expected one METADATA file in {wheel_path.relative_to(REPO_ROOT)}, found {len(metadata_paths)}"
            )
        metadata_text = zf.read(metadata_paths[0]).decode("utf-8", "replace")
    return Parser(policy=policy.default).parsestr(metadata_text)


def _normalize_package_name(name: str) -> str:
    return name.replace("_", "-").lower()


def _clean_metadata_value(value: str | None) -> str:
    if not value:
        return ""
    return " ".join(value.split())


def _source_url_from_metadata(metadata: Message) -> str:
    project_urls = metadata.get_all("Project-URL") or []
    preferred_labels = ("homepage", "repository", "source", "source code", "code", "documentation")
    parsed_urls: list[tuple[str, str]] = []
    for raw in project_urls:
        if "," not in raw:
            continue
        label, url = raw.split(",", 1)
        parsed_urls.append((label.strip().lower(), url.strip()))
    for preferred in preferred_labels:
        for label, url in parsed_urls:
            if label == preferred and url:
                return url
    if parsed_urls:
        return parsed_urls[0][1]
    return _clean_metadata_value(metadata.get("Home-page"))


def _license_from_metadata(package_name: str, metadata: Message) -> str:
    normalized_name = _normalize_package_name(package_name)
    if normalized_name in WHEEL_LICENSE_OVERRIDES:
        return WHEEL_LICENSE_OVERRIDES[normalized_name]
    license_expression = _clean_metadata_value(metadata.get("License-Expression"))
    if license_expression:
        return license_expression
    license_field = _clean_metadata_value(metadata.get("License"))
    if license_field and len(license_field) <= 120:
        return license_field
    for classifier in metadata.get_all("Classifier") or []:
        if classifier.startswith("License ::"):
            return classifier.rsplit("::", 1)[-1].strip()
    return license_field or "Not declared in wheel metadata"


def collect_wheel_notices() -> tuple[WheelNotice, ...]:
    notices: list[WheelNotice] = []
    for relative_path in PYODIDE_BUNDLED_WHEELS:
        wheel_path = WEB_ROOT / relative_path
        if not wheel_path.exists():
            raise FileNotFoundError(f"Missing bundled Python wheel: {wheel_path.relative_to(REPO_ROOT)}")
        metadata = _read_wheel_metadata(wheel_path)
        package_name = _clean_metadata_value(metadata.get("Name"))
        version = _clean_metadata_value(metadata.get("Version"))
        if not package_name or not version:
            raise RuntimeError(f"Wheel metadata is missing Name or Version: {wheel_path.relative_to(REPO_ROOT)}")
        notices.append(
            WheelNotice(
                filename=wheel_path.name,
                package_name=package_name,
                version=version,
                license_expression=_license_from_metadata(package_name, metadata),
                source_url=_source_url_from_metadata(metadata),
                bundled_path=relative_path.as_posix(),
            )
        )
    return tuple(
        sorted(notices, key=lambda notice: (_normalize_package_name(notice.package_name), notice.version))
    )


def _render_component_summary_rows(components: tuple[ComponentNotice, ...]) -> str:
    rows = []
    for component in components:
        rows.append(
            "                    <tr>\n"
            f"                        <td><strong>{_escape(component.display_name)}</strong></td>\n"
            f"                        <td>{_escape(component.version)}</td>\n"
            f"                        <td>{_license_html(component.license_expression, component.license_text_anchors)}</td>\n"
            f"                        <td><code>{_escape(component.bundled_path)}</code></td>\n"
            f"                        <td>{_external_link(component.source_url)}</td>\n"
            "                    </tr>"
        )
    return "\n".join(rows)


def _render_component_notices(components: tuple[ComponentNotice, ...]) -> str:
    rendered = []
    for component in components:
        component_id = f"component-{_slug(component.display_name)}"
        rendered.append(
            f'            <article class="component" id="{_escape(component_id)}">\n'
            f"                <h3>{_escape(component.display_name)}</h3>\n"
            "                <dl>\n"
            "                    <dt>Version</dt>\n"
            f"                    <dd>{_escape(component.version)}</dd>\n"
            "                    <dt>License</dt>\n"
            f"                    <dd>{_license_html(component.license_expression, component.license_text_anchors)}</dd>\n"
            "                    <dt>Bundled path</dt>\n"
            f"                    <dd><code>{_escape(component.bundled_path)}</code></dd>\n"
            "                    <dt>Notice</dt>\n"
            f"                    <dd>{_escape(component.notice)}</dd>\n"
            "                    <dt>Source</dt>\n"
            f"                    <dd>{_external_link(component.source_url)}</dd>\n"
            "                </dl>\n"
            "            </article>"
        )
    return "\n\n".join(rendered)


def _render_wheel_rows(wheels: tuple[WheelNotice, ...]) -> str:
    rows = []
    for wheel in wheels:
        anchors = _license_anchors_for_expression(wheel.license_expression)
        rows.append(
            "                    <tr>\n"
            f"                        <td><code>{_escape(wheel.bundled_path)}</code></td>\n"
            f"                        <td>{_escape(wheel.package_name)}</td>\n"
            f"                        <td>{_escape(wheel.version)}</td>\n"
            f"                        <td>{_license_html(wheel.license_expression, anchors)}</td>\n"
            f"                        <td>{_external_link(wheel.source_url)}</td>\n"
            "                    </tr>"
        )
    return "\n".join(rows)


def _render_license_sections(
    components: tuple[ComponentNotice, ...],
    wheels: tuple[WheelNotice, ...],
) -> str:
    anchors: set[str] = set()
    for component in components:
        anchors.update(component.license_text_anchors)
    for wheel in wheels:
        anchors.update(_license_anchors_for_expression(wheel.license_expression))

    ordered = [
        "license-mit",
        "license-apache-2.0",
        "license-mpl-2.0",
        "license-ofl-1.1",
        "license-bsd-3-clause",
    ]
    rendered = []
    for anchor in ordered:
        if anchor not in anchors:
            continue
        title, text_factory = LICENSE_TEXTS[anchor]
        rendered.append(
            f'            <article class="component license-anchor" id="{_escape(anchor)}">\n'
            f"                <h3>{_escape(title)}</h3>\n"
            f"                <pre>{_escape(_normalize_notice_text(text_factory()))}</pre>\n"
            "            </article>"
        )
    return "\n\n".join(rendered)


def render_open_source_notices() -> str:
    project_version = read_project_version()
    components = _component_notices(project_version)
    wheels = collect_wheel_notices()
    component_rows = _render_component_summary_rows(components)
    component_notices = _render_component_notices(components)
    wheel_rows = _render_wheel_rows(wheels)
    license_sections = _render_license_sections(components, wheels)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="Content-Security-Policy" content="
    default-src 'self';
    style-src 'self' 'unsafe-inline';
    font-src 'self' data:;
    img-src 'self' data:;
    object-src 'none';
    base-uri 'self';
    form-action 'self';
    ">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>gbdraw Open Source Notices</title>
    <style>
        @font-face {{
            font-family: 'Inter';
            src: url('./vendor/fonts/inter/inter-latin-400-normal.woff2') format('woff2');
            font-style: normal;
            font-weight: 400;
            font-display: swap;
        }}
        @font-face {{
            font-family: 'Inter';
            src: url('./vendor/fonts/inter/inter-latin-700-normal.woff2') format('woff2');
            font-style: normal;
            font-weight: 700;
            font-display: swap;
        }}
        @font-face {{
            font-family: 'Noto Sans JP';
            src: url('./vendor/fonts/noto-sans-jp/noto-sans-jp-japanese-400-normal.woff2') format('woff2');
            font-style: normal;
            font-weight: 400;
            font-display: swap;
        }}

        :root {{
            color-scheme: light;
            --bg: #f8fafc;
            --panel: #ffffff;
            --border: #cbd5e1;
            --text: #1e293b;
            --muted: #475569;
            --accent: #1d4ed8;
            --accent-soft: #dbeafe;
            --code: #0f172a;
        }}

        * {{ box-sizing: border-box; }}
        html {{ scroll-behavior: smooth; }}
        body {{
            margin: 0;
            background: linear-gradient(180deg, #f8fafc 0%, #eef2ff 100%);
            color: var(--text);
            font-family: 'Inter', 'Noto Sans JP', sans-serif;
            line-height: 1.6;
        }}

        a {{
            color: var(--accent);
            text-decoration: none;
        }}

        a:hover,
        a:focus-visible {{
            text-decoration: underline;
        }}

        .shell {{
            max-width: 1120px;
            margin: 0 auto;
            padding: 2rem 1rem 4rem;
        }}

        .hero,
        .component {{
            background: var(--panel);
            border: 1px solid rgba(148, 163, 184, 0.34);
            border-radius: 0.5rem;
            box-shadow: 0 14px 40px rgba(15, 23, 42, 0.06);
        }}

        .hero {{
            padding: 1.5rem;
        }}

        .eyebrow {{
            display: inline-block;
            font-size: 0.78rem;
            font-weight: 700;
            letter-spacing: 0.08em;
            text-transform: uppercase;
            background: var(--accent-soft);
            color: var(--accent);
            border-radius: 999px;
            padding: 0.3rem 0.7rem;
            margin-bottom: 0.9rem;
        }}

        h1,
        h2,
        h3 {{
            line-height: 1.2;
            margin: 0 0 0.8rem;
        }}

        h1 {{ font-size: clamp(2rem, 4vw, 3rem); }}
        h2 {{ font-size: 1.4rem; margin-top: 2rem; }}
        h3 {{ font-size: 1.05rem; margin-top: 1.3rem; }}

        p,
        li,
        td,
        th {{
            color: var(--muted);
        }}

        .actions {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.75rem;
            margin-top: 1.25rem;
        }}

        .button {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-height: 2.7rem;
            padding: 0.7rem 1rem;
            border-radius: 999px;
            border: 1px solid #93c5fd;
            background: #eff6ff;
            color: #1d4ed8;
            font-weight: 700;
        }}

        .button.secondary {{
            background: #ffffff;
            border-color: var(--border);
            color: var(--text);
        }}

        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 1rem;
            background: var(--panel);
            border-radius: 0.5rem;
            overflow: hidden;
        }}

        .summary-table th,
        .summary-table td {{
            text-align: left;
            padding: 0.85rem 0.9rem;
            border-bottom: 1px solid rgba(203, 213, 225, 0.75);
            vertical-align: top;
        }}

        .summary-table th {{
            background: #e2e8f0;
            color: #0f172a;
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}

        .summary-table tr:last-child td {{
            border-bottom: 0;
        }}

        .component {{
            margin-top: 1rem;
            padding: 1.1rem 1.15rem;
        }}

        .component dl {{
            display: grid;
            grid-template-columns: minmax(7rem, 10rem) 1fr;
            gap: 0.35rem 0.9rem;
            margin: 0;
        }}

        .component dt {{
            font-weight: 700;
            color: var(--text);
        }}

        .component dd {{
            margin: 0;
            color: var(--muted);
        }}

        code {{
            font-family: ui-monospace, SFMono-Regular, Consolas, 'Liberation Mono', monospace;
            background: #e2e8f0;
            color: var(--code);
            border-radius: 0.4rem;
            padding: 0.08rem 0.35rem;
        }}

        pre {{
            margin: 0.9rem 0 0;
            padding: 1rem;
            overflow-x: auto;
            white-space: pre-wrap;
            word-break: break-word;
            border-radius: 0.5rem;
            border: 1px solid rgba(148, 163, 184, 0.34);
            background: #0f172a;
            color: #e2e8f0;
            font-size: 0.83rem;
            line-height: 1.55;
            font-family: ui-monospace, SFMono-Regular, Consolas, 'Liberation Mono', monospace;
        }}

        .license-anchor {{
            scroll-margin-top: 2rem;
        }}

        .muted-note {{
            font-size: 0.92rem;
        }}
    </style>
</head>
<body>
    <main class="shell">
        <section class="hero">
            <div class="eyebrow">Open Source Notices</div>
            <h1>gbdraw Web App Licensing and Attribution</h1>
            <p>
                Licenses and attribution for gbdraw, its bundled browser assets, fonts, Pyodide runtime,
                and Python wheels used by the web app.
            </p>
            <p>
                Core runtime assets are served from local bundled files, so the app does not fetch Pyodide,
                Python wheels, icons, or fonts from third-party CDNs while it runs. The hosted gbdraw.app
                deployment separately loads the Google Analytics tag for aggregate page-usage metrics.
            </p>
            <div class="actions">
                <a class="button" href="./index.html">Back to gbdraw</a>
                <a class="button secondary" href="https://github.com/satoshikawato/gbdraw" target="_blank" rel="noopener noreferrer">Project Repository</a>
                <a class="button secondary" href="#license-mit">Jump to License Texts</a>
            </div>
        </section>

        <section>
            <h2>Bundled Components</h2>
            <table class="summary-table">
                <thead>
                    <tr>
                        <th scope="col">Component</th>
                        <th scope="col">Version / Provenance</th>
                        <th scope="col">License</th>
                        <th scope="col">Bundled Path</th>
                        <th scope="col">Source URL</th>
                    </tr>
                </thead>
                <tbody>
{component_rows}
                </tbody>
            </table>
            <p class="muted-note">
                The Tailwind Play CDN script is fetched by the vendor refresh tool from an upstream endpoint that does
                not expose an immutable version string in the downloaded bundle. gbdraw vendors that exact snapshot
                locally and serves the local copy at runtime.
            </p>
        </section>

        <section>
            <h2>Component Notices</h2>
{component_notices}
        </section>

        <section>
            <h2>Bundled Python and Pyodide Wheels</h2>
            <table class="summary-table">
                <thead>
                    <tr>
                        <th scope="col">Bundled Wheel</th>
                        <th scope="col">Package</th>
                        <th scope="col">Version</th>
                        <th scope="col">License</th>
                        <th scope="col">Source URL</th>
                    </tr>
                </thead>
                <tbody>
{wheel_rows}
                </tbody>
            </table>
            <p class="muted-note">
                Package-specific license files, when present upstream, remain embedded in the wheel archives listed
                above. License expressions shown here prefer explicit local overrides for ambiguous metadata, then
                <code>License-Expression</code>, then <code>License</code>, then license classifiers.
            </p>
        </section>

        <section>
            <h2>License Texts</h2>
{license_sections}
        </section>
    </main>
</body>
</html>
"""


def write_open_source_notices() -> None:
    NOTICES_PATH.write_text(render_open_source_notices(), encoding="utf-8")


def check_open_source_notices() -> bool:
    expected = render_open_source_notices()
    current = NOTICES_PATH.read_text(encoding="utf-8") if NOTICES_PATH.exists() else ""
    if current == expected:
        return True
    diff = difflib.unified_diff(
        current.splitlines(keepends=True),
        expected.splitlines(keepends=True),
        fromfile=str(NOTICES_PATH.relative_to(REPO_ROOT)),
        tofile="generated open-source notices",
    )
    sys.stderr.writelines(diff)
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate gbdraw web open source notices from local metadata.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit non-zero if gbdraw/web/open-source-notices.html is stale.",
    )
    args = parser.parse_args(argv)
    if args.check:
        return 0 if check_open_source_notices() else 1
    write_open_source_notices()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
