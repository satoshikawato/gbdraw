"""The bcbio-gff guards.

gbdraw reads GFF3 with bcbio-gff, and the Web app runs one vendored wheel of it.
``test_bcbio_gff_version_is_the_same_everywhere`` tells a maintainer which files to
update together when that version changes. ``test_mixed_separators_need_gbdraws_repair``
tells a maintainer when the installed bcbio-gff has started to split mixed ``" ; "`` and
``"; "`` attribute separators itself, and when gbdraw's own repair has stopped working.
"""

from __future__ import annotations

import io
import json
import re
import sys
from pathlib import Path

from BCBio import GFF
import pytest

from gbdraw._web_assets import PYODIDE_LOCAL_WHEELS
from gbdraw.io.input_text import gff3_stream, text_lines

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parents[1]
MIXED_SEPARATORS_GFF3 = REPO_ROOT / "tests" / "fixtures" / "input_variants" / "ggn93_exonerate_style.gff3"


def _web_wheel() -> tuple[str, str]:
    """The path and the version of the one bcbio-gff wheel the Web app runs."""

    wheels = [path for path in PYODIDE_LOCAL_WHEELS if path.name.startswith("bcbio_gff-")]
    assert len(wheels) == 1, f"PYODIDE_LOCAL_WHEELS must hold exactly one bcbio-gff wheel: {wheels}"
    # PEP 427 file name: {distribution}-{version}-{python tag}-{abi tag}-{platform tag}.whl
    return wheels[0].as_posix(), wheels[0].name.split("-")[1]


def test_bcbio_gff_version_is_the_same_everywhere() -> None:
    wheel_path, version = _web_wheel()

    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    requirements = [
        requirement
        for requirement in pyproject["project"]["dependencies"]
        if re.match(r"bcbio-gff\b", requirement)
    ]
    assert requirements == [f"bcbio-gff=={version}"], (
        f"pyproject.toml must require bcbio-gff=={version}, the version of the Web wheel "
        f"{wheel_path} in gbdraw/_web_assets.py; update pyproject.toml, recipe/meta.yaml, "
        "gbdraw/_web_assets.py and tools/immutable-web-assets.json together"
    )

    recipe = (REPO_ROOT / "recipe" / "meta.yaml").read_text(encoding="utf-8")
    assert re.search(rf"^\s+- bcbio-gff =={re.escape(version)}\s*$", recipe, re.MULTILINE), (
        f"recipe/meta.yaml must list the run requirement '- bcbio-gff =={version}', the version "
        "of the Web wheel; update it together with pyproject.toml"
    )

    immutable = json.loads((REPO_ROOT / "tools" / "immutable-web-assets.json").read_text(encoding="utf-8"))
    assert wheel_path in immutable["files"], (
        f"tools/immutable-web-assets.json must list {wheel_path} under 'files'; update it "
        "together with gbdraw/_web_assets.py"
    )


def test_mixed_separators_need_gbdraws_repair() -> None:
    """bcbio-gff stops on mixed separators; gbdraw's repair lets the file parse.

    When the first assertion fails, the installed bcbio-gff splits mixed separators
    itself and ``_gff3_attributes`` in gbdraw/io/input_text.py may be removed. When
    the second fails, gbdraw's repair is broken.
    """

    text = MIXED_SEPARATORS_GFF3.read_text(encoding="utf-8")
    with pytest.raises(AssertionError):
        list(GFF.parse(io.StringIO(text)))

    with MIXED_SEPARATORS_GFF3.open(encoding="utf-8") as handle:
        records = list(GFF.parse(gff3_stream(text_lines(handle))))
    assert records
    assert any(record.features for record in records)
