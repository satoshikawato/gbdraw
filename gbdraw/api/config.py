"""Config helpers for the public API layer."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from typing import Mapping

from gbdraw.config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.config.toml import load_config_toml  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]


def load_default_config() -> dict:
    """Load the packaged default config.toml as a dict."""

    return load_config_toml("gbdraw.data", "config.toml")


def apply_config_overrides(
    config: GbdrawConfig | dict | None,
    overrides: Mapping[str, object] | None,
) -> GbdrawConfig:
    """Apply CLI-style overrides and return a typed GbdrawConfig.

    Notes:
        - If `config` is None, the packaged default config.toml is loaded.
        - If `config` is a GbdrawConfig, it is converted to a dict (modeled fields only).
    """

    if config is None:
        config_dict = load_default_config()
    elif isinstance(config, GbdrawConfig):
        # asdict only includes modeled fields; this is intended for typed configs.
        config_dict = asdict(config)
    else:
        config_dict = deepcopy(config)

    if overrides:
        try:
            config_dict = modify_config_dict(config_dict, **overrides)
        except TypeError as exc:
            raise ValidationError(f"Invalid config override: {exc}") from exc

    return GbdrawConfig.from_dict(config_dict)


__all__ = ["apply_config_overrides", "load_default_config", "GbdrawConfig"]
