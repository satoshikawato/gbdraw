"""Versioned capabilities exposed by the Python runtime to the Web application.

The browser must negotiate this manifest once when its diagram worker starts.
Capabilities are derived from their owning Python constants; callers must not
infer support by inspecting implementation source or command-line help text.
"""

from __future__ import annotations

from typing import Any

from gbdraw.features.shapes import FEATURE_RENDERING_VALUES
from gbdraw.session_io import (
    DEPTH_FILE_ENCODING,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
)
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
    UNKNOWN_FIELD_POLICY,
)
from gbdraw.tracks import (
    SUPPORTED_CIRCULAR_TRACK_RENDERERS,
    SUPPORTED_LINEAR_TRACK_RENDERERS,
)


WEB_RUNTIME_CAPABILITY_SCHEMA = 1
WEB_RENDER_PROTOCOL = 2
WEB_RENDER_OPTIONS_SCHEMA = 1


def get_web_runtime_capabilities() -> dict[str, Any]:
    """Return the JSON-compatible contract implemented by this runtime.

    A fresh object is returned so a caller cannot mutate process-wide capability
    state. Increment ``WEB_RUNTIME_CAPABILITY_SCHEMA`` for a structural manifest
    change and ``WEB_RENDER_PROTOCOL`` for an incompatible worker/render contract.
    """

    return {
        "schema": WEB_RUNTIME_CAPABILITY_SCHEMA,
        "renderProtocol": WEB_RENDER_PROTOCOL,
        "request": {
            "currentSchema": CANONICAL_REQUEST_SCHEMA,
            "supportedSchemas": sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS),
            "unknownFieldPolicy": UNKNOWN_FIELD_POLICY,
        },
        "resources": {
            "encodings": ["base64", DEPTH_FILE_ENCODING],
        },
        "rendering": {
            "optionSchema": WEB_RENDER_OPTIONS_SCHEMA,
            "featureRenderings": sorted(FEATURE_RENDERING_VALUES),
            "circularTrackRenderers": sorted(
                SUPPORTED_CIRCULAR_TRACK_RENDERERS
            ),
            "linearTrackRenderers": sorted(SUPPORTED_LINEAR_TRACK_RENDERERS),
        },
        "analysis": {
            "proteinLosatCacheSchema": PROTEIN_LOSAT_CACHE_SCHEMA,
            "nucleotideLosatCacheSchema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
            "losatDerivedCacheSchema": LOSAT_DERIVED_CACHE_SCHEMA,
            "proteinIdentityManifestSchema": PROTEIN_IDENTITY_MANIFEST_SCHEMA,
        },
    }


__all__ = [
    "WEB_RENDER_OPTIONS_SCHEMA",
    "WEB_RENDER_PROTOCOL",
    "WEB_RUNTIME_CAPABILITY_SCHEMA",
    "get_web_runtime_capabilities",
]
