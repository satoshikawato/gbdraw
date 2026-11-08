#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import logging
import sys
from functools import lru_cache
from pathlib import Path

# tomllib is available in Python 3.11+; use tomli as fallback for 3.10
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from importlib import resources

# Traversable lives in importlib.resources.abc since 3.11; Python 3.14 removed importlib.abc.Traversable
if sys.version_info >= (3, 11):
    from importlib.resources.abc import Traversable
else:
    from importlib.abc import Traversable

from typing import cast

from ..exceptions import ConfigError

logger = logging.getLogger(__name__)


@lru_cache(maxsize=8)
def _config_text(path: Path) -> str:
    """The text of one config file, read and logged once per process: the packaged file does not change."""
    logger.info("Loading config file: %s", path)
    return path.read_bytes().decode("utf-8")


def load_config_toml(config_directory: str, config_file: str) -> dict:
    absolute_config_path = None  # Initialize outside of try for scope in exception block
    try:
        # Generate the path object for the 'config.toml' file
        config_path: Traversable = resources.files(config_directory).joinpath(config_file)
        # Convert the path to an absolute path
        # Installed packages are directories, so files() returns a Path here.
        absolute_config_path = cast("Path", config_path).resolve()
        # Parsed on every call, so each caller gets a dict of its own to change.
        config_dict = tomllib.loads(_config_text(absolute_config_path))
    except FileNotFoundError as e:
        raise ConfigError(
            f"Configuration file missing: {absolute_config_path}. "
            "Ensure gbdraw is installed with package data."
        ) from e
    except Exception as e:
        raise ConfigError(
            f"Failed to load configuration from {absolute_config_path}: {e}"
        ) from e
    return config_dict


__all__ = ["load_config_toml"]


