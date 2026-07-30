"""Support helpers shared by the browser UI and offline web asset tooling."""

from .capabilities import (
    WEB_RENDER_OPTIONS_SCHEMA,
    WEB_RENDER_PROTOCOL,
    WEB_RUNTIME_CAPABILITY_SCHEMA,
    get_web_runtime_capabilities,
)
__all__ = [
    "WEB_RENDER_OPTIONS_SCHEMA",
    "WEB_RENDER_PROTOCOL",
    "WEB_RUNTIME_CAPABILITY_SCHEMA",
    "get_web_runtime_capabilities",
    "render_canonical_web_request",
    "render_embedded_canonical_web_request",
]


def __getattr__(name: str) -> object:
    if name in {
        "render_canonical_web_request",
        "render_embedded_canonical_web_request",
    }:
        from . import request_render

        return getattr(request_render, name)
    raise AttributeError(name)
