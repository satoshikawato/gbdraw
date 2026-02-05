#!/usr/bin/env python
# coding: utf-8

import logging
import sys

# tomllib is available in Python 3.11+; use tomli as fallback for 3.10
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from importlib import resources
from importlib.abc import Traversable

from ..exceptions import ConfigError

logger = logging.getLogger(__name__)


def load_config_toml(config_directory: str, config_file: str) -> dict:
    absolute_config_path = None  # Initialize outside of try for scope in exception block
    try:
        # Generate the path object for the 'config.toml' file
        config_path: Traversable = resources.files(config_directory).joinpath(config_file)
        # Convert the path to an absolute path
        absolute_config_path = config_path.resolve()  # type: ignore
        # Display or log the absolute path
        logger.info(f"INFO: Loading config file: {absolute_config_path}")
        # Open the file and load the configuration
        with open(absolute_config_path, "rb") as config_toml:
            config_dict = tomllib.load(config_toml)
    except FileNotFoundError as e:
        logger.error(f"Failed to load configs from {absolute_config_path}: {e}")
        raise ConfigError(
            f"Configuration file missing: {absolute_config_path}. "
            "Ensure gbdraw is installed with package data."
        ) from e
    except Exception as e:
        logger.error(f"Failed to load configs from {absolute_config_path}: {e}")
        raise ConfigError(
            f"Failed to load configuration from {absolute_config_path}: {e}"
        ) from e
    return config_dict


__all__ = ["load_config_toml"]


