"""SE-08(a): derived label-filtering maps are never preserved as user settings."""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path

import pandas as pd

from gbdraw.api.config import load_default_config
from gbdraw.api.diagram import _resolve_diagram_options_config
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.labels.filtering import (
    DERIVED_LABEL_FILTERING_KEYS,
    get_label_text,
    preprocess_label_filtering,
)
from gbdraw.web_support.config_overrides import validate_and_project_web_config_overrides

REPO_ROOT = Path(__file__).resolve().parents[1]
SESSIONS = REPO_ROOT / "tests" / "fixtures" / "sessions"
PROVENANCE = json.loads(
    (SESSIONS / "se08-main-linear-cli.provenance.json").read_text(encoding="utf-8")
)


def _main_cli_sidecar_config() -> dict:
    raw = gzip.decompress(
        (SESSIONS / "se08-main-linear-cli.v42.gbdraw-session.json.gz").read_bytes()
    )
    session = json.loads(raw)
    return session["renderRequest"]["diagramOptions"]["config"]


def test_derived_key_set_matches_the_preprocessor_outputs():
    filtering = preprocess_label_filtering({"blacklist_keywords": []})
    assert DERIVED_LABEL_FILTERING_KEYS == frozenset(filtering) - {"blacklist_keywords"}


def test_main_cli_linear_session_preserves_no_raw_filtering():
    config = _main_cli_sidecar_config()
    assert DERIVED_LABEL_FILTERING_KEYS <= set(config["labels"]["filtering"])
    projected = validate_and_project_web_config_overrides(
        mode="linear",
        config=config,
        overrides={},
        managed_paths=PROVENANCE["managedPaths"],
    )
    assert projected == {}


def test_main_web_saved_raw_with_derived_maps_is_dropped_on_load_and_save():
    stored = PROVENANCE["webStoredUnmanagedConfigOverrides"]
    assert DERIVED_LABEL_FILTERING_KEYS <= set(stored["labels.filtering.raw"])
    projected = validate_and_project_web_config_overrides(
        mode="linear",
        overrides=stored,
        managed_paths=PROVENANCE["managedPaths"],
        require_unmanaged_only=True,
    )
    assert projected == {}


def test_tobacco_gallery_full_config_preserves_no_raw_filtering():
    session = json.loads(
        (REPO_ROOT / "gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json")
        .read_text(encoding="utf-8")
    )
    config = session["renderRequest"]["diagramOptions"]["config"]
    projected = validate_and_project_web_config_overrides(
        mode="circular",
        config=config,
        managed_paths=["labels.filtering.blacklist_keywords"],
    )
    assert "labels.filtering.raw" not in projected


def test_real_unmanaged_raw_change_is_kept_without_derived_maps():
    raw = {
        "blacklist_keywords": [],
        "qualifier_priority": {"gene": ["gene"]},
        "whitelist_map": None,
        "priority_map": {"CDS": ["product"]},
        "label_override_rules": None,
    }
    projected = validate_and_project_web_config_overrides(
        mode="linear",
        overrides={"labels.filtering.raw": raw},
        managed_paths=["labels.filtering.blacklist_keywords"],
        require_unmanaged_only=True,
    )
    assert projected["labels.filtering.raw"] == {
        "blacklist_keywords": [],
        "qualifier_priority": {"gene": ["gene"]},
    }
    config = load_default_config()
    config["labels"]["filtering"] = copy.deepcopy(raw)
    projected = validate_and_project_web_config_overrides(
        mode="linear",
        config=config,
        managed_paths=["labels.filtering.blacklist_keywords"],
    )
    assert projected["labels.filtering.raw"] == {
        "blacklist_keywords": [],
        "qualifier_priority": {"gene": ["gene"]},
    }


class _Feature:
    type = "CDS"
    qualifiers = {"product": ["NADH"], "locus_tag": ["R0"]}


def test_attached_tables_replace_stale_compiled_maps():
    config = _main_cli_sidecar_config()
    config["labels"]["filtering"]["priority_map"] = {"CDS": ["product"]}
    config["labels"]["filtering"]["whitelist_map"] = {"CDS": {"product": ["PRIVATE"]}}
    for options_type in (LinearDiagramOptions, CircularDiagramOptions):
        options = options_type(
            config=copy.deepcopy(config),
            qualifier_priority_table=pd.DataFrame(
                [{"feature_type": "CDS", "priorities": "locus_tag"}]
            ),
        )
        filtering = preprocess_label_filtering(
            _resolve_diagram_options_config(options).labels.filtering.as_dict()
        )
        assert filtering["priority_map"] == {"CDS": ["locus_tag"]}
        assert get_label_text(_Feature(), filtering) == "R0"
