"""Real GC content, GC skew and depth configurators with chosen fields, for tests of the series builders."""

from __future__ import annotations

from typing import TypeVar

from gbdraw.api.config import apply_config_overrides
from gbdraw.config.models import CircularRenderProfile
from gbdraw.configurators import DepthConfigurator, GcContentConfigurator, GcSkewConfigurator
from gbdraw.io.colors import load_default_colors

_Config = TypeVar("_Config", GcContentConfigurator, GcSkewConfigurator, DepthConfigurator)


def _with(config: _Config, fields: dict[str, object]) -> _Config:
    for name, value in fields.items():
        assert hasattr(config, name), name
        setattr(config, name, value)
    return config


def _profile() -> CircularRenderProfile:
    return CircularRenderProfile(apply_config_overrides(None, {}))


def gc_content_config(**fields: object) -> GcContentConfigurator:
    config = GcContentConfigurator(
        window=100, step=10, dinucleotide="GC", profile=_profile(), default_colors_df=load_default_colors("", "default")
    )
    return _with(config, fields)


def gc_skew_config(**fields: object) -> GcSkewConfigurator:
    config = GcSkewConfigurator(
        window=100, step=10, dinucleotide="GC", profile=_profile(), default_colors_df=load_default_colors("", "default")
    )
    return _with(config, fields)


def depth_config(**fields: object) -> DepthConfigurator:
    return _with(DepthConfigurator(100, 10, _profile()), fields)
