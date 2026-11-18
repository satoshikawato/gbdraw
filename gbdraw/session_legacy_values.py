"""Values that older Sessions and schema 1-2 requests wrote, and what they read as now.

The canonical-request decoder (:mod:`gbdraw.session_request_codec`), the
Session migration of a draft (:mod:`gbdraw.session_io`) and the replay of a
recorded command line (:mod:`gbdraw.cli_utils.legacy_session_args`) read the
same renames from here, and the values of an option in a recorded command line
(``_option_value``, ``_option_all_values``).
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Mapping, Sequence

# The feature types an older writer drew when it named none.
LEGACY_FEATURE_TYPES: tuple[str, ...] = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "repeat_region",
)
# The shapes an older writer drew these feature types as when it left the shape unset.
LEGACY_FEATURE_SHAPES: Mapping[str, str] = MappingProxyType({"repeat_region": "rectangle"})
# The Circular multi-record size modes an older writer wrote under another name.
LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES: Mapping[str, str] = MappingProxyType({"sqrt": "auto"})
# The depth and GC content axes of a config written before they had a label style (v0.14.5 and earlier) drew a
# number on each large tick: the config sections, by name, their override paths, and the style they read as
# (D-16, D-24, 2026-11-08). A stored full config and stored config overrides without the style read so; every
# writer of config overrides since states it.
LEGACY_SCALAR_AXIS_LABEL_SECTIONS: tuple[str, ...] = ("depth", "gc_content")
SCALAR_AXIS_LABEL_PATHS: tuple[str, ...] = tuple(
    f"objects.{section}.axis_labels" for section in LEGACY_SCALAR_AXIS_LABEL_SECTIONS
)
LEGACY_SCALAR_AXIS_LABELS = "numbered"
# The resource IDs under which a CLI or Python render of gbdraw (request schemas 2 and 5-9) stored a comparison
# table it was given (`--comparisons_table`, a `LinearComparison`): this prefix and its number. Those writers stored
# such a table as a protein table; the Web named its LOSATP tables otherwise (OV-3573).
LEGACY_GIVEN_COMPARISON_RESOURCE_PREFIX = "comparison-explicit-"
# The request schemas of those writers: 2 and 5-9 on main from 2026-07-18, 9 in the releases.
LEGACY_GIVEN_COMPARISON_SCHEMAS: tuple[int, ...] = (2, 5, 6, 7, 8, 9)
# The column every protein table those writers stored has. Under those schemas a table under that prefix and without
# this column holds no protein matches; every other stored protein table holds them. Python's decoder and the Web Load
# read them by this rule (OV-3512, OV-3573, OV-3574).
LEGACY_GIVEN_COMPARISON_PROTEIN_COLUMN = "query_protein_id"


def _option_all_values(args: Sequence[str], *names: str) -> list[str]:
    values: list[str] = []
    for index, token in enumerate(args):
        text = str(token)
        for name in names:
            if text == name:
                if index + 1 < len(args):
                    values.append(str(args[index + 1]))
                break
            prefix = f"{name}="
            if text.startswith(prefix):
                values.append(text[len(prefix):])
                break
    return values


def _option_value(args: Sequence[str], *names: str) -> str | None:
    for index, token in enumerate(args):
        text = str(token)
        for name in names:
            prefix = f"{name}="
            if text.startswith(prefix):
                return text[len(prefix):]
            if text == name and index + 1 < len(args):
                return str(args[index + 1])
    return None


__all__ = [
    "LEGACY_GIVEN_COMPARISON_PROTEIN_COLUMN",
    "LEGACY_GIVEN_COMPARISON_RESOURCE_PREFIX",
    "LEGACY_GIVEN_COMPARISON_SCHEMAS",
    "LEGACY_SCALAR_AXIS_LABELS",
    "LEGACY_SCALAR_AXIS_LABEL_SECTIONS",
    "SCALAR_AXIS_LABEL_PATHS",
    "LEGACY_FEATURE_SHAPES",
    "LEGACY_FEATURE_TYPES",
    "LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES",
]
