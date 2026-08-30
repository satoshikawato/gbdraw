"""Schema-owned validation for Web-preserved configuration leaves."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any, Literal

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.config.modify import (
    canonical_config_override_paths,
    config_to_raw_dict,
    validate_config_overrides,
)
from gbdraw.config.models import GbdrawConfig
from gbdraw.exceptions import ValidationError


def _raw_config_leaf(config: Mapping[str, Any], path: str) -> object:
    if path == "labels.filtering.raw":
        return config["labels"]["filtering"]
    current: object = config
    for key in path.split("."):
        if not isinstance(current, Mapping) or key not in current:
            raise KeyError(path)
        current = current[key]
    return current


def _validated_base_config(config: object) -> GbdrawConfig:
    if config is None:
        return apply_config_overrides(None, None)
    if not isinstance(config, Mapping):
        raise ValidationError("diagramOptions.config must be an object or null.")
    try:
        return GbdrawConfig.from_dict(config)
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ValidationError("diagramOptions.config is invalid.") from exc


def validate_and_project_web_config_overrides(
    *,
    mode: Literal["circular", "linear"],
    config: object = None,
    overrides: object = None,
    managed_paths: object = None,
    require_unmanaged_only: bool = False,
) -> dict[str, object]:
    """Return validated effective leaves that the active Web GUI does not own.

    The typed configuration model remains the validity owner. ``managed_paths``
    describes only the current GUI projection; it cannot make a path valid.
    """

    if mode not in {"circular", "linear"}:
        raise ValidationError(f"Unsupported Web configuration mode: {mode!r}.")
    if overrides is not None and not isinstance(overrides, Mapping):
        raise ValidationError("diagramOptions.configOverrides must be an object or null.")
    if managed_paths is None:
        managed_paths = []
    if not isinstance(managed_paths, list) or not all(
        isinstance(path, str) for path in managed_paths
    ):
        raise ValidationError("Managed config override paths must be a string array.")

    canonical_paths = canonical_config_override_paths()
    unknown_managed = sorted(set(managed_paths) - canonical_paths)
    if unknown_managed:
        raise ValidationError(
            "The Web GUI declares unknown managed config path(s): "
            + ", ".join(unknown_managed)
            + "."
        )
    managed = frozenset(managed_paths)
    explicit = dict(overrides or {})
    validated = dict(validate_config_overrides(explicit))

    if require_unmanaged_only:
        dual_owned = sorted(set(validated) & managed)
        if dual_owned:
            raise ValidationError(
                "Preserved config override path is GUI-managed: "
                + ", ".join(dual_owned)
                + "."
            )

    options_type = CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    options_type(config_overrides=validated)
    base = _validated_base_config(config)
    try:
        effective = apply_config_overrides(base, validated)
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ValidationError("Config override values are incompatible.") from exc

    default_raw = config_to_raw_dict(apply_config_overrides(None, None))
    effective_raw = config_to_raw_dict(effective)
    candidates = {
        path
        for path in canonical_paths
        if _raw_config_leaf(effective_raw, path) != _raw_config_leaf(default_raw, path)
    }
    candidates.update(validated)

    projected: dict[str, object] = {}
    for path in sorted(candidates - managed):
        if path in validated:
            projected[path] = validated[path]
        else:
            value = _raw_config_leaf(effective_raw, path)
            try:
                options_type(config_overrides={path: value})
            except ValidationError:
                # Full configs contain both mode branches. Only derived leaves
                # accepted by the active typed options become preserved overrides.
                continue
            projected[path] = value

    # Re-run mode validation against values originating from a full config.
    options_type(config_overrides=projected)
    return projected


def validate_web_config_overrides_json(
    mode: object,
    config_json: object,
    overrides_json: object,
    managed_paths_json: object,
    require_unmanaged_only: object = False,
) -> str:
    """JSON boundary used by the existing diagram-Worker helper channel."""

    try:
        config = json.loads(str(config_json))
        overrides = json.loads(str(overrides_json))
        managed_paths = json.loads(str(managed_paths_json))
        projected = validate_and_project_web_config_overrides(
            mode=str(mode),  # type: ignore[arg-type]
            config=config,
            overrides=overrides,
            managed_paths=managed_paths,
            require_unmanaged_only=bool(require_unmanaged_only),
        )
    except ValidationError as exc:
        return json.dumps({"error": str(exc)}, ensure_ascii=False)
    except (KeyError, TypeError, ValueError):
        return json.dumps(
            {"error": "Web config override payload is invalid."},
            ensure_ascii=False,
        )
    return json.dumps({"overrides": projected}, ensure_ascii=False)


__all__ = [
    "validate_and_project_web_config_overrides",
    "validate_web_config_overrides_json",
]
