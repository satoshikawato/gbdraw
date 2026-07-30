"""Config helpers for the public API layer."""

from __future__ import annotations

from typing import Mapping

from gbdraw.config.modify import (  # type: ignore[reportMissingImports]
    _apply_validated_config_overrides,
    config_to_raw_dict,
)
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.config.toml import load_config_toml  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError


def load_default_config() -> dict:
    """Load the packaged default config.toml as a dict."""

    return load_config_toml("gbdraw.data", "config.toml")


def apply_config_overrides(
    config: GbdrawConfig | dict | None,
    overrides: Mapping[str, object] | None,
) -> GbdrawConfig:
    """Apply canonical dotted-path overrides and return a typed GbdrawConfig.

    Notes:
        - If `config` is None, the packaged default config.toml is loaded.
        - If `config` is a GbdrawConfig, it is converted to a dict (modeled fields only).
    """

    if config is None:
        config_dict = load_default_config()
    elif isinstance(config, GbdrawConfig):
        if not overrides:
            return config
        config_dict = config_to_raw_dict(config)
    else:
        try:
            base_config = GbdrawConfig.from_dict(config)
        except ValidationError:
            raise
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError(f"Invalid configuration: {exc}") from exc
        if not overrides:
            return base_config
        config_dict = config_to_raw_dict(base_config)

    if overrides:
        config_dict = _apply_validated_config_overrides(config_dict, overrides)

    return GbdrawConfig.from_dict(config_dict)


__all__ = ["apply_config_overrides", "load_default_config", "GbdrawConfig"]
